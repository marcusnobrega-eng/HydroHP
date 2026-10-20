function metrics = test_hydrostatic_reconstruction
%TEST_HYDROSTATIC_RECONSTRUCTION P0.5: positivity and balance of the
% first-order hydrostatic reconstruction (Audusse et al., 2004).

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
options = {'include_friction', false, 'cfl_target', 0.8, ...
    'reconstruction', "hydrostatic"};

% Lake at rest over a raised bed and a contraction, and around a bump whose
% crest stands above the water surface (dry cells inside the domain).
cases = {[10; 6; 10], [0; 0.4; 0], 1.5; [10; 10; 10], [0; 2.0; 0], 1.0};
dry_cell_count = zeros(1, 2);
for index = 1:size(cases, 1)
    [cells, faces, dx_m] = make_grid(cases{index, 1}, cases{index, 2}, 40);
    geometry = hydrohp_compile_geometry(cells, faces);
    area = areas_at_stage(cells, cases{index, 3});
    initial_area = area;
    discharge = zeros(size(area));
    adjustment_m3 = 0;
    for step = 1:300
        [area, discharge, diagnostic] = hydrohp_varying_step_compiled( ...
            geometry, area, discharge, dx_m, 0.02, options{:});
        adjustment_m3 = adjustment_m3 + diagnostic.volume_adjustment_m3;
    end
    assert(max(abs(area - initial_area)) < 1e-10);
    assert(max(abs(discharge)) < 1e-10);
    assert(adjustment_m3 == 0);
    dry_cell_count(index) = sum(initial_area == 0);
end
assert(dry_cell_count(1) == 0 && dry_cell_count(2) > 0);

% On a flat bed both reconstructions coincide exactly.
[cells, faces, dx_m] = make_grid([8; 8; 8], [0; 0; 0], 60);
geometry = hydrohp_compile_geometry(cells, faces);
x = [cells.x_m].';
area = areas_at_stage(cells, 1.0 * (x < 40));
discharge = zeros(size(area));
area_stage = area;
discharge_stage = discharge;
for step = 1:100
    [area, discharge] = hydrohp_varying_step_compiled(geometry, area, ...
        discharge, dx_m, 0.01, options{:});
    [area_stage, discharge_stage] = hydrohp_varying_step_compiled(geometry, ...
        area_stage, discharge_stage, dx_m, 0.01, ...
        include_friction=false, cfl_target=0.8);
end
assert(isequal(area, area_stage) && isequal(discharge, discharge_stage));

% The P0.1 steep, widening channel creates volume with stage reconstruction
% but none with hydrostatic reconstruction, and no raw area turns negative.
capture_root = string(tempname);
mkdir(capture_root);
cleanup = onCleanup(@() rmdir(capture_root, 's'));
input = steep_channel_input(fullfile(capture_root, "capture.mat"));
input.reconstruction = "hydrostatic";
steep = hydrohp_run_varying(input);
assert(~isfile(input.negative_area_capture_file));
assert(steep.initialization.ledger.wet_dry_adjustment_m3 == 0);
assert(steep.ledger.wet_dry_adjustment_m3 == 0);
assert(abs(steep.initialization.ledger.relative_residual) < 1e-12);
assert(steep.metadata.reconstruction == "hydrostatic");

% Uniform flow on a mild slope converges at first order and is more
% accurate than stage reconstruction at the same spacing.
error_hydrostatic = zeros(1, 2);
cell_counts = [20, 40];
for index = 1:2
    error_hydrostatic(index) = uniform_depth_error(cell_counts(index), ...
        "hydrostatic");
end
error_stage = uniform_depth_error(20, "stage");
assert(error_hydrostatic(1) < error_stage);
assert(error_hydrostatic(2) < 0.6 * error_hydrostatic(1));

% The reference kernel has no hydrostatic option.
input.kernel = "reference";
input.negative_area_capture_file = "";
try
    hydrohp_run_varying(input);
    error('hydrohp:TestFailed', 'Hydrostatic reconstruction ran on the reference kernel.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:InvalidVaryingInput'));
end

metrics = struct( ...
    'emerging_bump_dry_cells', dry_cell_count(2), ...
    'uniform_depth_error_hydrostatic_dx30', error_hydrostatic(1), ...
    'uniform_depth_error_hydrostatic_dx15', error_hydrostatic(2), ...
    'uniform_depth_error_stage_dx30', error_stage);
fprintf(['HydroHP hydrostatic-reconstruction checks passed: uniform-flow ' ...
    'depth error %.1f%% (dx 30 m) and %.1f%% (dx 15 m), stage %.1f%%.\n'], ...
    100 * error_hydrostatic(1), 100 * error_hydrostatic(2), 100 * error_stage);
end

function relative_error = uniform_depth_error(cell_count, reconstruction)
slope = 0.003;
sections = channel_sections([8; 8; 8], slope);
input = struct( ...
    'section_data', sections, ...
    'face_x_m', linspace(0, 600, cell_count + 1).', ...
    'hydrograph_time_s', [0; 1e5], ...
    'hydrograph_discharge_m3s', [2; 2], ...
    'duration_s', 1, ...
    'downstream', struct('type', "normal_depth", 'slope', slope), ...
    'base_discharge_m3s', 2, ...
    'initialization_type', "normal_depth", ...
    'spinup_duration_s', 1800, ...
    'steady_tolerance_per_s', 1e-4, ...
    'maximum_dt_s', 5, ...
    'cfl_target', 0.7, ...
    'reconstruction', reconstruction);
result = hydrohp_run_varying(input);
cells = hydrohp_map_sections_to_grid(sections, result.x_m);
normal_area = interp1(cells(1).table.conveyance_m3s, cells(1).table.area_m2, ...
    2 / sqrt(slope), 'pchip');
normal_depth = hydrohp_depth_from_area(cells(1), normal_area);
inner = find(result.x_m > 100 & result.x_m < 300);
depth = arrayfun(@(k) hydrohp_depth_from_area(cells(k), ...
    result.initialization.final_area_m2(k)), inner);
relative_error = max(abs(depth - normal_depth)) / normal_depth;
end

function input = steep_channel_input(capture_file)
input = struct( ...
    'section_data', channel_sections([8; 8; 40], 0.03), ...
    'face_x_m', linspace(0, 600, 21).', ...
    'hydrograph_time_s', [0; 600], ...
    'hydrograph_discharge_m3s', [2; 2], ...
    'duration_s', 60, ...
    'downstream', struct('type', "normal_depth", 'slope', 0.002), ...
    'base_discharge_m3s', 2, ...
    'initialization_type', "normal_depth", ...
    'spinup_duration_s', 600, ...
    'steady_tolerance_per_s', 1e9, ...
    'maximum_dt_s', 5, ...
    'cfl_target', 0.7, ...
    'negative_area_capture_file', capture_file);
end

function sections = channel_sections(bottom_width_m, slope)
chainage_m = [0; 300; 600];
rows = table();
for index = 1:numel(chainage_m)
    w = bottom_width_m(index);
    invert = 20 - slope * chainage_m(index);
    section = table( ...
        repmat("S" + index, 4, 1), repmat(chainage_m(index), 4, 1), ...
        [0; 2; 2 + w; 4 + w], invert + [4; 0; 0; 4], ...
        0.035 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section]; %#ok<AGROW>
end
sections = hydrohp_build_section_tables(rows, 0.02, "scm");
end

function [cells, faces, dx_m] = make_grid(width_m, invert_m, cell_count)
chainage_m = [0; 50; 100];
rows = table();
for section_index = 1:numel(chainage_m)
    section_rows = table( ...
        repmat("S" + section_index, 4, 1), ...
        repmat(chainage_m(section_index), 4, 1), ...
        width_m(section_index) * [0; 0.1; 0.9; 1], ...
        invert_m(section_index) + [3; 0; 0; 3], ...
        0.03 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section_rows]; %#ok<AGROW>
end
surveyed = hydrohp_build_section_tables(rows, 0.02, "scm");
face_x = linspace(0, 100, cell_count + 1).';
dx_m = face_x(2) - face_x(1);
cell_x = face_x(1:end-1) + dx_m / 2;
cells = hydrohp_map_sections_to_grid(surveyed, cell_x);
faces = hydrohp_map_sections_to_grid(surveyed, face_x);
end

function area = areas_at_stage(sections, stage_m)
stage_m = stage_m .* ones(numel(sections), 1);
area = zeros(numel(sections), 1);
for index = 1:numel(sections)
    if stage_m(index) > sections(index).invert_elevation_m
        area(index) = interp1( ...
            sections(index).table.water_surface_elevation_m, ...
            sections(index).table.area_m2, stage_m(index), 'pchip');
    end
end
end
