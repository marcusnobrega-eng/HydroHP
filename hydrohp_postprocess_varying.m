function files = hydrohp_postprocess_varying(result, input, options)
%HYDROHP_POSTPROCESS_VARYING V5 post-processing products for a varying-section run.
%
% files = hydrohp_postprocess_varying(result, input) writes the figures,
% videos, and CSV tables of the V5 post-processing (post_processing.m and the
% scripts it calls) for a run of HYDROHP_RUN_VARYING. File names, layouts,
% colours, and formatting follow the V5 scripts. The function reads the result
% and the run input and changes neither. files lists the paths it wrote.
%
% Differences from V5, all forced by the varying-section data:
% - The kernel does not store its time step, so the top width replaces the
%   Courant number in the summary charts, the states video, and the CSVs.
% - Values are at cell centres. Each cross-section panel shows the surveyed
%   section nearest to the inlet, middle, or outlet cell, filled to the water
%   level interpolated at its chainage.
% - Observation points are the run probes when every requested point is a
%   probe (every accepted step), otherwise cells of the stored history.

arguments
    result (1,1) struct
    input (1,1) struct
    options.output_folder (1,1) string = "Modeling_Results"
    options.run_id (1,1) string = "1"
    options.run_name (1,1) string = "HydroHP"
    options.observation_x_m (:,1) double = zeros(0, 1)
    options.observation_labels (:,1) string = strings(0, 1)
    options.figures (1,1) logical = true
    options.videos (1,1) logical = true
    options.csv (1,1) logical = true
    options.csv_layout (1,1) string {mustBeMember(options.csv_layout, ...
        ["by_section", "by_state"])} = "by_section"
    options.arrival_depth_m (1,1) double {mustBeFinite, mustBeNonnegative} = 0.5
    options.person_height_m (1,1) double {mustBeFinite, mustBePositive} = 1.70
    options.person_width_m (1,1) double {mustBeFinite, mustBePositive} = 0.35
    options.specific_weight_n_m3 (1,1) double {mustBeFinite, mustBePositive} = 9810
    options.computational_time_s (1,1) double = NaN
    options.frame_rate (1,1) double {mustBeFinite, mustBePositive} = 20
end

if ~isfield(result, 'time_s') || numel(result.time_s) < 2 || ...
        ~isfield(result, 'area_m2') || ~isfield(input, 'section_data')
    error('hydrohp:NoHistory', ...
        'Post-processing needs a result with at least two stored times and its run input.');
end
g = 9.81;
if isfield(input, 'gravity_m_s2')
    g = input.gravity_m_s2;
end
minimum_area = 1e-10;
if isfield(input, 'minimum_area_m2')
    minimum_area = input.minimum_area_m2;
end

x = result.x_m(:);
time_s = result.time_s(:);
cells = hydrohp_map_sections_to_grid(input.section_data, x);
bed = [cells.invert_elevation_m].';
states = cell_states(result.area_m2, result.discharge_m3s, cells, g, minimum_area);
force = person_force(states.depth_m, states.velocity_m_s, options, g);
observation = observation_series(result, cells, options, g, minimum_area);

folder = options.output_folder;
if ~isfolder(folder)
    mkdir(folder);
end
prefix = options.run_id + "_" + options.run_name + "_";
files = struct();

if options.figures
    files.surf_plots = write_surf_plots(folder, prefix, x, time_s, states);
    files.summary_charts = write_summary_charts(folder, prefix, x, bed, ...
        cells, states, observation);
    files.maximum_profiles = write_maximum_profiles(folder, prefix, x, bed, ...
        result.maxima, max(force, [], 1).', observation.x_m);
end
if options.videos
    files.dam_breach_video = write_dam_breach_video(folder, prefix, x, bed, ...
        time_s, states, force, observation.x_m, options.frame_rate);
    files.states_video = write_states_video(folder, prefix, x, time_s, ...
        states, options.frame_rate);
    files.cross_section_video = write_cross_section_video(folder, prefix, x, ...
        bed, time_s, states, cells, input.section_data, options.frame_rate);
end
if options.csv
    files.detailed_output = write_detailed_output(folder, prefix, x, time_s, ...
        states, options.csv_layout);
    files.detailed_output_xs = write_detailed_output_xs(folder, prefix, x, ...
        time_s, states);
    files.dam_breach_outputs = write_dam_breach_outputs(folder, prefix, x, ...
        input, observation, options);
end
end

% ---------------------------------------------------------------------------
% States

function states = cell_states(area, discharge, cells, g, minimum_area)
% Rows are times, columns are cells. Depth follows the solver (pchip of the
% area table) and top width is linear in depth, as in the kernel.
depth = zeros(size(area));
width = zeros(size(area));
for index = 1:numel(cells)
    table_data = cells(index).table;
    local_area = min(max(area(:, index), 0), table_data.area_m2(end));
    depth(:, index) = interp1(table_data.area_m2, table_data.depth_m, ...
        local_area, 'pchip');
    width(:, index) = interp1(table_data.depth_m, table_data.top_width_m, ...
        depth(:, index), 'linear');
end
wet = area > minimum_area;
velocity = zeros(size(area));
velocity(wet) = discharge(wet) ./ area(wet);
froude = zeros(size(area));
wet = wet & width > 0;
froude(wet) = abs(velocity(wet)) ./ sqrt(g * area(wet) ./ width(wet));
states = struct('area_m2', area, 'discharge_m3s', discharge, ...
    'depth_m', depth, 'velocity_m_s', velocity, 'froude', froude, ...
    'top_width_m', width, 'stage_m', depth + [cells.invert_elevation_m]);
end

function force_tf = person_force(depth, velocity, options, g)
% Same formula and units as post_processing_Dam_Breach.m, so the V5 and
% varying-section products agree.
gamma = options.specific_weight_n_m3;
width = options.person_width_m;
height = options.person_height_m;
area_person = width * max(depth, height);
velocity_force = (gamma / g) * velocity.^2 / (2 * g) .* area_person / 1000 / 10;
pressure = gamma * depth.^2 / 2 - gamma * max(depth - height, 0).^2 / 2;
force_tf = velocity_force + pressure * width / 1000 / 10;
end

function observation = observation_series(result, cells, options, g, minimum_area)
x = result.x_m(:);
has_probes = isfield(result, 'probes') && isfield(result.probes, 'x_m') && ...
    ~isempty(result.probes.x_m);
observation_x = options.observation_x_m;
if isempty(observation_x)
    if has_probes
        observation_x = result.probes.x_m(:);
    else
        observation_x = x([1, ceil(numel(x) / 2), numel(x)]);
    end
end
if numel(observation_x) > 5
    error('hydrohp:TooManyObservationPoints', ...
        'The V5 products hold at most five observation points.');
end
labels = options.observation_labels;
if isempty(labels)
    labels = compose("x = %g m", observation_x);
elseif numel(labels) ~= numel(observation_x)
    error('hydrohp:InvalidObservationLabels', ...
        'Give one observation label per observation point.');
end

probe_index = zeros(numel(observation_x), 1);
if has_probes
    for index = 1:numel(observation_x)
        match = find(abs(result.probes.x_m(:) - observation_x(index)) <= ...
            1e-9 * max(1, abs(observation_x(index))), 1, 'first');
        if ~isempty(match)
            probe_index(index) = match;
        end
    end
end
if all(probe_index > 0)
    cell_index = result.probes.cell_index(probe_index);
    time_s = result.probes.time_s(:);
    area = result.probes.area_m2(:, probe_index);
    discharge = result.probes.discharge_m3s(:, probe_index);
    source = "probes";
else
    cell_index = zeros(numel(observation_x), 1);
    for index = 1:numel(observation_x)
        [~, cell_index(index)] = min(abs(x - observation_x(index)));
    end
    time_s = result.time_s(:);
    area = result.area_m2(:, cell_index);
    discharge = result.discharge_m3s(:, cell_index);
    source = "history";
end
observation = cell_states(area, discharge, cells(cell_index), g, minimum_area);
observation.force_tf = person_force(observation.depth_m, ...
    observation.velocity_m_s, options, g);
observation.time_s = time_s;
observation.x_m = observation_x;
observation.cell_index = cell_index(:);
observation.labels = labels;
observation.source = source;
end

% ---------------------------------------------------------------------------
% Figures (post_processing.m and post_processing_Dam_Breach.m)

function path = write_surf_plots(folder, prefix, x, time_s, states)
fig = new_figure([2, 0, 8, 10]);
cleanup = onCleanup(@() close(fig));
time_h = time_s / 3600;
panels = {states.froude, 'Froude Number'; states.depth_m, 'y (m)'; ...
    states.stage_m, 'WSE (m)'};
for index = 1:3
    ax = subplot(3, 1, index, 'Parent', fig);
    surf(ax, x, time_h, panels{index, 1});
    view(ax, 0, 90);
    color_bar = colorbar(ax);
    colormap(ax, 'jet');
    color_bar.TickDirection = 'out';
    shading(ax, 'interp');
    xlabel(ax, 'Distance from the dam [m]', 'Interpreter', 'latex');
    ylabel(ax, 'Elapsed time [h]', 'Interpreter', 'latex');
    ylabel(color_bar, panels{index, 2}, 'Interpreter', 'latex');
    xlim(ax, [x(1), x(end)]);
    ylim(ax, [0, time_h(end)]);
    set(ax, 'FontName', 'Garamond', 'FontSize', 12, 'LineWidth', 1.5, ...
        'TickLength', [0.02 0.01], 'TickDir', 'out');
    box(ax, 'on');
end
path = fullfile(folder, prefix + "Surf_Plots.pdf");
exportgraphics(fig, path, 'ContentType', 'image', 'Colorspace', 'rgb', ...
    'Resolution', 600);
end

function path = write_summary_charts(folder, prefix, x, bed, cells, states, observation)
fig = new_figure([2, 0, 8, 10]);
cleanup = onCleanup(@() close(fig));
line_style = {'-', '--', '-.', ':', '-'};
line_width = [2 1.75 1.5 1 .75];
time_min = observation.time_s / 60;
panels = {observation.discharge_m3s, 'Flow Discharge (m\textsuperscript{3}/s)'; ...
    observation.velocity_m_s, 'Velocity (m/s)'; ...
    observation.depth_m, 'Water Depths (m)'; ...
    observation.froude, 'Froude Number'; ...
    observation.top_width_m, 'Top Width (m)'};
for index = 1:5
    ax = subplot(3, 2, index, 'Parent', fig);
    hold(ax, 'on');
    for point = 1:numel(observation.labels)
        plot(ax, time_min, panels{index, 1}(:, point), ...
            'LineStyle', line_style{point}, 'LineWidth', line_width(point), ...
            'Color', 'k');
    end
    xlabel(ax, 'Elapsed time (min)', 'Interpreter', 'latex');
    ylabel(ax, panels{index, 2}, 'Interpreter', 'latex');
    legend(ax, observation.labels, 'Interpreter', 'latex', 'Location', 'best');
end

% Rating curve: inlet and middle cells against the Manning normal depth of
% the inlet cell, as in V5.
middle = ceil(numel(x) / 2);
dx_m = x(2) - x(1);
slope = (bed(1) - bed(2)) / dx_m;
if slope <= 0
    slope = (bed(1) - bed(end)) / (x(end) - x(1));
end
inlet_table = cells(1).table;
maximum_depth = max(states.depth_m(:));
normal_depth = linspace(0, min(1.1 * maximum_depth, inlet_table.depth_m(end)), 200).';
normal_discharge = NaN(size(normal_depth));
if slope > 0
    normal_discharge = interp1(inlet_table.depth_m, inlet_table.conveyance_m3s, ...
        normal_depth, 'pchip') * sqrt(slope);
end
ax = subplot(3, 2, 6, 'Parent', fig);
hold(ax, 'on');
plot(ax, states.discharge_m3s(2:end, 1), states.depth_m(2:end, 1), ...
    'LineStyle', '--', 'LineWidth', 2, 'Color', 'k');
plot(ax, states.discharge_m3s(2:end, middle), states.depth_m(2:end, middle), ...
    'LineStyle', ':', 'LineWidth', 2, 'Color', 'k');
plot(ax, normal_discharge, normal_depth, 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'k');
xlabel(ax, 'Flow Discharge (m\textsuperscript{3}/s)', 'Interpreter', 'latex');
ylabel(ax, 'Water Depth (m)', 'Interpreter', 'latex');
ylim(ax, [0, upper_limit(1.1 * maximum_depth)]);
legend(ax, 'Q(Inlet)', 'Q(Nx/2)', '$Q_{n}$ (L)', 'Interpreter', 'latex', ...
    'Location', 'best');
path = fullfile(folder, prefix + "Summary_Charts.pdf");
% Probe series hold every step, so MATLAB suggests an image; V5 keeps vector.
warning_state = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
restore = onCleanup(@() warning(warning_state));
exportgraphics(fig, path, 'ContentType', 'vector');
end

function path = write_maximum_profiles(folder, prefix, x, bed, maxima, ...
        maximum_force, observation_x)
fig = new_figure([0, 0, 7, 12]);
cleanup = onCleanup(@() close(fig));
colors = legacy_colors();
x_km = x / 1000;
maximum_stage = maxima.stage_m(:);
maximum_velocity = maxima.velocity_m_s(:);
observation_bed = interp1(x, bed, observation_x, 'linear', 'extrap');

ax = subplot(3, 1, 1, 'Parent', fig);
hold(ax, 'on');
plot(ax, x_km, bed, 'LineWidth', 4, 'LineStyle', '-', 'Color', 'black');
plot(ax, x_km, maximum_stage, 'LineWidth', 2, 'LineStyle', '-', ...
    'Color', colors.water);
fill(ax, [x_km; flipud(x_km)], [bed; flipud(maximum_stage)], colors.water);
xlabel(ax, 'Distance from the dam [km]', 'Interpreter', 'latex');
ylabel(ax, 'Water surface Elevation [m]', 'Interpreter', 'latex');
ylim(ax, stage_limits(bed, maximum_stage));
grid(ax, 'on');
title(ax, 'Maximum Water Surface Elevation', 'Interpreter', 'latex', 'FontSize', 12);
observation_lines(ax, observation_x, observation_bed, max(maximum_stage));

ax = subplot(3, 1, 2, 'Parent', fig);
hold(ax, 'on');
plot(ax, x_km, maximum_velocity, 'LineWidth', 2, 'LineStyle', '-', ...
    'Color', colors.velocity);
xlabel(ax, 'Distance from the dam [km]', 'Interpreter', 'latex');
ylabel(ax, 'Velocity [$\mathrm{m~s^{-1}}$]', 'Interpreter', 'latex');
ylim(ax, [0, upper_limit(max(maximum_velocity))]);
grid(ax, 'on');
title(ax, 'Maximum Velocity', 'Interpreter', 'latex', 'FontSize', 12);
observation_lines(ax, observation_x, zeros(size(observation_x)), ...
    max(maximum_velocity));

ax = subplot(3, 1, 3, 'Parent', fig);
hold(ax, 'on');
plot(ax, x_km, maximum_force, 'LineWidth', 2, 'LineStyle', '-', ...
    'Color', colors.force);
xlabel(ax, 'Distance from the dam [km]', 'Interpreter', 'latex');
ylabel(ax, 'Force [tf]', 'Interpreter', 'latex');
ylim(ax, [0, upper_limit(max(maximum_force))]);
grid(ax, 'on');
title(ax, 'Force [tf]', 'Interpreter', 'latex', 'FontSize', 12);
observation_lines(ax, observation_x, zeros(size(observation_x)), ...
    max(maximum_force));
path = fullfile(folder, prefix + "Maximum_WSE_Velocity_Depth.pdf");
exportgraphics(fig, path, 'ContentType', 'vector');
end

% ---------------------------------------------------------------------------
% Videos (post_processing_Dam_Breach.m, states_post_processing.m, and
% cross_section_post_processing.m)

function path = write_dam_breach_video(folder, prefix, x, bed, time_s, ...
        states, force, observation_x, frame_rate)
path = fullfile(folder, prefix + "Dam_Breach_Like_.mp4");
[writer, fig, cleanup] = open_video(path, [0, 0, 7, 12], frame_rate); %#ok<ASGLU>
colors = legacy_colors();
x_km = x / 1000;
zero_ground = zeros(size(x));
stage_range = stage_limits(bed, states.stage_m(:));
maximum_velocity = upper_limit(max(states.velocity_m_s(:)));
maximum_force = upper_limit(max(force(:)));
observation_bed = interp1(x, bed, observation_x, 'linear', 'extrap');
for frame = 1:numel(time_s)
    clf(fig);
    sgtitle(fig, sprintf('Time = %.2f (min)', time_s(frame) / 60), ...
        'FontSize', 18, 'Interpreter', 'latex');
    stage = states.stage_m(frame, :).';
    velocity = states.velocity_m_s(frame, :).';
    frame_force = force(frame, :).';

    ax = subplot(3, 1, 1, 'Parent', fig);
    hold(ax, 'on');
    plot(ax, x_km, bed, 'LineWidth', 4, 'LineStyle', '-', 'Color', 'black');
    plot(ax, x_km, stage, 'LineWidth', 2, 'LineStyle', '-', 'Color', colors.water);
    fill(ax, [x_km; flipud(x_km)], [bed; flipud(stage)], colors.water);
    xlabel(ax, 'Distance from the dam [km]', 'Interpreter', 'latex');
    ylabel(ax, 'Water surface Elevation [m]', 'Interpreter', 'latex');
    ylim(ax, stage_range);
    grid(ax, 'on');
    title(ax, 'Water Surface Elevation [m]', 'Interpreter', 'latex');
    observation_lines(ax, observation_x, observation_bed, stage_range(2));

    ax = subplot(3, 1, 2, 'Parent', fig);
    hold(ax, 'on');
    plot(ax, x_km, velocity, 'LineWidth', 2, 'LineStyle', '-', ...
        'Color', colors.velocity);
    fill(ax, [x_km; flipud(x_km)], [zero_ground; flipud(velocity)], colors.velocity);
    xlabel(ax, 'Distance from the dam [km]', 'Interpreter', 'latex');
    ylabel(ax, 'Velocity [m/s]', 'Interpreter', 'latex');
    ylim(ax, [0, maximum_velocity]);
    grid(ax, 'on');
    title(ax, 'Velocity [m/s]', 'Interpreter', 'latex');
    observation_lines(ax, observation_x, zeros(size(observation_x)), maximum_velocity);

    ax = subplot(3, 1, 3, 'Parent', fig);
    hold(ax, 'on');
    plot(ax, x_km, frame_force, 'LineWidth', 2, 'LineStyle', '-', ...
        'Color', colors.force);
    fill(ax, [x_km; flipud(x_km)], [zero_ground; flipud(frame_force)], colors.force);
    xlabel(ax, 'Distance from the dam [km]', 'Interpreter', 'latex');
    ylabel(ax, 'Force [tf]', 'Interpreter', 'latex');
    ylim(ax, [0, maximum_force]);
    grid(ax, 'on');
    title(ax, 'Force [tf]', 'Interpreter', 'latex');
    observation_lines(ax, observation_x, zeros(size(observation_x)), maximum_force);
    writeVideo(writer, getframe(fig));
end
end

function path = write_states_video(folder, prefix, x, time_s, states, frame_rate)
path = fullfile(folder, prefix + "Model_Results.mp4");
[writer, fig, cleanup] = open_video(path, [2, 0.5, 12, 7.5], frame_rate); %#ok<ASGLU>
hec_ras_colors = [52/255 85/255 132/255; 0 1 1; 0 128/255 1; 0 255/255 0; ...
    1 1 0; 1 128/255 0; 1 0 0; 128/255 0 128/255];
panels = {states.depth_m, '$h(t)~\mathrm{[m]}$', 'Depth'; ...
    states.discharge_m3s, '$q_1(t)~\mathrm{[m^3/s]}$', 'Discharge'; ...
    states.area_m2, '$q_2(t)~\mathrm{[m^2]}$', 'Flow Area'; ...
    states.froude, '$F_r(t)~\mathrm{[-]}$', 'Froude'; ...
    states.top_width_m, '$B(t)~\mathrm{[m]}$', 'Top Width'; ...
    states.velocity_m_s, '$v(t)~\mathrm{[m/s]}$', 'Velocity'};
limits = zeros(6, 2);
for index = 1:6
    values = panels{index, 1}(:);
    limits(index, :) = [min(values), max(values)];
    if limits(index, 2) <= limits(index, 1)
        limits(index, 2) = limits(index, 1) + 1;
    end
end
for frame = 1:numel(time_s)
    clf(fig);
    sgtitle(fig, sprintf('Time = %.0f (sec)', time_s(frame)), ...
        'FontSize', 18, 'Interpreter', 'latex');
    for index = 1:6
        ax = subplot(2, 3, index, 'Parent', fig);
        hold(ax, 'on');
        values = panels{index, 1}(frame, :).';
        plot(ax, x, values, 'k-', 'LineWidth', 2);
        fill(ax, [x; flipud(x)], [zeros(size(x)); flipud(values)], ...
            hec_ras_colors(index, :));
        [peak, peak_index] = max(values);
        scatter(ax, x(peak_index), peak, 50, 'MarkerFaceColor', ...
            hec_ras_colors(index, :), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        xlabel(ax, 'Distance (m)', 'Interpreter', 'latex');
        ylabel(ax, panels{index, 2}, 'Interpreter', 'latex');
        title(ax, panels{index, 3});
        ylim(ax, limits(index, :));
        grid(ax, 'on');
        set(ax, 'FontName', 'Garamond', 'FontSize', 12, 'FontWeight', 'Bold', ...
            'LineWidth', 1.5, 'TickLength', [0.02 0.01], 'TickDir', 'out');
    end
    writeVideo(writer, getframe(fig));
end
end

function path = write_cross_section_video(folder, prefix, x, bed, time_s, ...
        states, cells, section_data, frame_rate)
path = fullfile(folder, prefix + "Depth_WSE_Top_Width.mp4");
[writer, fig, cleanup] = open_video(path, [2, 0, 8, 10], frame_rate); %#ok<ASGLU>
colors = legacy_colors();
tick_size = [0.015 0.01];
chainage = [section_data.chainage_m].';
targets = x([1, ceil(numel(x) / 2), numel(x)]);
shown = zeros(3, 1);
for index = 1:3
    [~, shown(index)] = min(abs(chainage - targets(index)));
end
section_width = max(arrayfun(@(item) item.table.top_width_m(end), cells));
stage_range = stage_limits(bed, states.stage_m(:));
for frame = 1:numel(time_s)
    clf(fig);
    sgtitle(fig, sprintf('Time = %.0f (sec)', time_s(frame)), ...
        'FontSize', 18, 'Interpreter', 'latex');
    stage = states.stage_m(frame, :).';
    for index = 1:3
        section = section_data(shown(index));
        query = min(max(section.chainage_m, x(1)), x(end));
        section_stage = interp1(x, stage, query, 'linear');
        depth = max(section_stage - section.invert_elevation_m, 0);
        ax = subplot(3, 3, index, 'Parent', fig);
        hold(ax, 'on');
        plot(ax, section.points.offset_m, section.points.elevation_m, '-k', ...
            'LineWidth', 2, 'Marker', '*');
        fill_water(ax, section.points.offset_m, section.points.elevation_m, ...
            section_stage, colors.water);
        box(ax, 'on');
        if index == 1
            ylabel(ax, 'Elevation [m]', 'Interpreter', 'latex');
        end
        xlabel(ax, 'Station [m]', 'Interpreter', 'latex');
        title(ax, sprintf('x = %0.2f m, h = %0.2f m', section.chainage_m, depth), ...
            'FontSize', 12, 'Interpreter', 'latex');
        set(ax, 'FontSize', 12, 'FontName', 'Garamond', ...
            'TickLength', tick_size, 'TickDir', 'in');
    end

    width = states.top_width_m(frame, :).';
    right_margin = section_width / 2 + width / 2;
    left_margin = section_width / 2 - width / 2;
    ax = subplot(3, 3, [4 5 6], 'Parent', fig);
    hold(ax, 'on');
    plot(ax, x, right_margin, 'k', 'LineWidth', 2);
    plot(ax, x, left_margin, 'k', 'LineWidth', 2);
    fill(ax, [x; flipud(x)], [left_margin; flipud(right_margin)], colors.water);
    set(ax, 'YDir', 'reverse');
    xlabel(ax, '$x$ [m]', 'Interpreter', 'latex');
    ylabel(ax, 'Station [m]', 'Interpreter', 'latex');
    ylim(ax, [0, upper_limit(section_width)]);
    grid(ax, 'on');
    title(ax, sprintf('$B_{{max}}(t)$ = %0.2f m', max(width)), ...
        'FontSize', 16, 'Interpreter', 'latex');
    set(ax, 'FontSize', 12, 'FontName', 'Garamond', ...
        'TickLength', tick_size, 'TickDir', 'in');

    ax = subplot(3, 3, [7 8 9], 'Parent', fig);
    hold(ax, 'on');
    plot(ax, x, bed, 'LineWidth', 4, 'LineStyle', '-', 'Color', 'k');
    plot(ax, x, stage, 'LineWidth', 2, 'LineStyle', '-', 'Color', colors.water);
    fill(ax, [x; flipud(x)], [bed; flipud(stage)], colors.water);
    xlabel(ax, '$x$ [m]', 'Interpreter', 'latex');
    ylabel(ax, 'Water Surface Elevation [m]', 'Interpreter', 'latex');
    ylim(ax, stage_range);
    grid(ax, 'on');
    title(ax, sprintf('$WSE_{{max}}(t)$ = %0.2f m', max(stage)), ...
        'FontSize', 16, 'Interpreter', 'latex');
    set(ax, 'FontSize', 12, 'FontName', 'Garamond', ...
        'TickLength', tick_size, 'TickDir', 'in');
    writeVideo(writer, getframe(fig));
end
end

% ---------------------------------------------------------------------------
% CSV tables (Detailed_Output_Script.m and Detailed_Output_Dam_Breach_Like.m)

function path = write_detailed_output(folder, prefix, x, time_s, states, layout)
names = ["Depth (m)", "Discharge (m^3/s)", "Velocity (m/s)", "Froude (-)", ...
    "Top Width (m)", "Water Surface Elevation (m)"];
data = state_stack(states);
[time_count, cell_count, state_count] = size(data);
if layout == "by_section"
    values = reshape(permute(data, [1 3 2]), time_count, []);
    header_names = repmat(names, 1, cell_count);
    header_x = repelem(x.', state_count);
else
    values = reshape(data, time_count, []);
    header_names = repelem(names, cell_count);
    header_x = repmat(x.', 1, state_count);
end
headers = compose("%s, x(m) = %0.2f", header_names(:), header_x(:)).';
output = array2table([time_s, round(values, 3)], ...
    'VariableNames', ["Time (sec)", headers]);
path = fullfile(folder, prefix + "Detailed_Output.csv");
writetable(output, path, 'Delimiter', ',');
end

function path = write_detailed_output_xs(folder, prefix, x, time_s, states)
% One block per cell, separated by an empty row, as in HEC-RAS tables.
data = state_stack(states);
[time_count, cell_count, state_count] = size(data);
values = NaN((time_count + 1) * cell_count, state_count + 2);
for index = 1:cell_count
    rows = (index - 1) * (time_count + 1) + (1:time_count);
    values(rows, :) = [repmat(x(index), time_count, 1), time_s, ...
        reshape(data(:, index, :), time_count, state_count)];
end
output = array2table(values, 'VariableNames', {'x(m)', 'Elapsed Time (sec)', ...
    'Depth (m)', 'Discharge (m3/s)', 'Velocity (m/s)', 'Froude (-)', ...
    'Top Width (m)', 'WSE (m)'});
path = fullfile(folder, prefix + "Detailed_Output_XS.csv");
writetable(output, path, 'Delimiter', ',');
end

function path = write_dam_breach_outputs(folder, prefix, x, input, observation, options)
count = 5;
time_max_depth = NaN(1, count);
maximum_depth = NaN(1, count);
time_max_velocity = NaN(1, count);
maximum_velocity = NaN(1, count);
arrival_time = NaN(1, count);
maximum_force = NaN(1, count);
time_min = observation.time_s / 60;
for point = 1:numel(observation.x_m)
    depth = observation.depth_m(:, point);
    velocity = observation.velocity_m_s(:, point);
    [maximum_depth(point), position] = max(depth);
    time_max_depth(point) = time_min(position);
    [maximum_velocity(point), position] = max(velocity);
    time_max_velocity(point) = time_min(position);
    position = find(depth > depth(1) + options.arrival_depth_m, 1, 'first');
    if ~isempty(position)
        arrival_time(point) = time_min(position);
    end
    maximum_force(point) = max(observation.force_tf(:, point));
end
slot = compose("[%d]", (1:count).');
headers = ["Computational Time (min)", "Simulation Time (min)", ...
    "Spatial Discretization (m)", ...
    "Arrival time of maximum depth (min) " + slot.', ...
    "Maximum Depth (m) " + slot.', ...
    "Arrival time of maximum velocity " + slot.', ...
    "Maximum velocity (m/s) " + slot.', ...
    "Flood arrival time (min) " + slot.', ...
    "Maximum force (tf) " + slot.'];
values = [options.computational_time_s / 60, input.duration_s / 60, x(2) - x(1), ...
    time_max_depth, maximum_depth, time_max_velocity, maximum_velocity, ...
    arrival_time, maximum_force];
output = array2table(values, 'VariableNames', headers);
path = fullfile(folder, prefix + "Dam_Breach_Like_Outputs.csv");
writetable(output, path, 'Delimiter', ',');
end

% ---------------------------------------------------------------------------
% Helpers

function data = state_stack(states)
data = cat(3, states.depth_m, states.discharge_m3s, states.velocity_m_s, ...
    states.froude, states.top_width_m, states.stage_m);
end

function fig = new_figure(position_in)
fig = figure('Visible', 'off', 'Units', 'inches', 'Position', position_in, ...
    'Color', 'w');
end

function [writer, fig, cleanup] = open_video(path, position_in, frame_rate)
writer = VideoWriter(char(path), 'MPEG-4');
writer.Quality = 100;
writer.FrameRate = frame_rate;
open(writer);
fig = new_figure(position_in);
cleanup = onCleanup(@() close_video(writer, fig));
end

function close_video(writer, fig)
close(writer);
if isvalid(fig)
    close(fig);
end
end

function colors = legacy_colors()
colors = struct('water', [21, 179, 196] / 255, ...
    'velocity', [255, 99, 71] / 256, 'force', [25, 25, 112] / 256);
end

function limits = stage_limits(bed, stage)
limits = [0.98 * min(bed), max(1.01 * stage)];
if limits(2) <= limits(1)
    limits(2) = limits(1) + 1;
end
end

function value = upper_limit(value)
if ~(value > 0)
    value = 1;
end
end

function observation_lines(ax, observation_x, bottom, top)
for index = 1:numel(observation_x)
    plot(ax, [observation_x(index), observation_x(index)] / 1000, ...
        [bottom(index), top], 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black');
end
end

function fill_water(ax, offset, elevation, stage, color)
% Fill each wet stretch of a profile up to the water level, with the
% crossings of the water level inserted between profile points.
offset = offset(:).';
elevation = elevation(:).';
x_dense = offset(1);
z_dense = elevation(1);
for index = 1:numel(offset) - 1
    z0 = elevation(index);
    z1 = elevation(index + 1);
    if (z0 - stage) * (z1 - stage) < 0
        x_dense(end + 1) = offset(index) + (stage - z0) * ...
            (offset(index + 1) - offset(index)) / (z1 - z0); %#ok<AGROW>
        z_dense(end + 1) = stage; %#ok<AGROW>
    end
    x_dense(end + 1) = offset(index + 1); %#ok<AGROW>
    z_dense(end + 1) = z1; %#ok<AGROW>
end
wet = z_dense <= stage;
edges = diff([false, wet, false]);
first = find(edges == 1);
last = find(edges == -1) - 1;
for run = 1:numel(first)
    span = first(run):last(run);
    if numel(span) > 1 && any(z_dense(span) < stage)
        fill(ax, [x_dense(span), fliplr(x_dense(span))], ...
            [z_dense(span), stage * ones(1, numel(span))], color);
    end
end
end
