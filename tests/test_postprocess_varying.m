function metrics = test_postprocess_varying
%TEST_POSTPROCESS_VARYING The V5 post-processing products for a varying-
% section run: every file is written, the CSV values equal the solver states,
% and observation points use probes when the run has them.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
input = struct( ...
    'section_data', channel([1; 0.5; 0]), ...
    'face_x_m', linspace(0, 600, 25).', ...
    'hydrograph_time_s', [0; 60; 300; 600; 1800], ...
    'hydrograph_discharge_m3s', [0; 30; 10; 0; 0], ...
    'duration_s', 1800, ...
    'downstream', struct('type', "normal_depth", 'slope', 1/600), ...
    'initialization_type', "dry", ...
    'maximum_dt_s', 2, ...
    'output_interval_s', 300, ...
    'probe_x_m', [0; 300; 600]);
result = hydrohp_run_varying(input);
folder = string(tempname);
cleanup = onCleanup(@() remove_folder(folder));

files = hydrohp_postprocess_varying(result, input, 'output_folder', folder, ...
    'run_id', "T", 'run_name', "Test");
expected = ["surf_plots", "summary_charts", "maximum_profiles", ...
    "dam_breach_video", "states_video", "cross_section_video", ...
    "detailed_output", "detailed_output_xs", "dam_breach_outputs"];
assert(isequal(sort(string(fieldnames(files))), sort(expected(:))));
for name = expected
    information = dir(files.(name));
    assert(isscalar(information) && information.bytes > 0, ...
        'Missing or empty output %s.', name);
    assert(startsWith(information.name, "T_Test_"));
end

% Detailed output: one row per stored time, six states per cell.
cells = hydrohp_map_sections_to_grid(input.section_data, result.x_m);
depth = zeros(size(result.area_m2));
for index = 1:numel(cells)
    depth(:, index) = hydrohp_depth_from_area(cells(index), result.area_m2(:, index));
end
detailed = readmatrix(files.detailed_output);
cell_count = numel(result.x_m);
assert(isequal(size(detailed), [numel(result.time_s), 1 + 6 * cell_count]));
assert(max(abs(detailed(:, 1) - result.time_s)) == 0);
assert(max(abs(detailed(:, 2:6:end) - round(depth, 3)), [], 'all') == 0);
assert(max(abs(detailed(:, 7:6:end) - ...
    round(depth + [cells.invert_elevation_m], 3)), [], 'all') == 0);
header = string(detectImportOptions(files.detailed_output).VariableNames);
assert(numel(header) == 1 + 6 * cell_count);

by_state = hydrohp_postprocess_varying(result, input, 'output_folder', folder, ...
    'run_id', "S", 'figures', false, 'videos', false, 'csv_layout', "by_state");
detailed_by_state = readmatrix(by_state.detailed_output);
assert(max(abs(detailed_by_state(:, 1 + (1:cell_count)) - round(depth, 3)), ...
    [], 'all') == 0);

% Cross-section table: one block per cell and an empty separator row.
xs = readmatrix(files.detailed_output_xs);
assert(size(xs, 2) == 8);
assert(size(xs, 1) >= (numel(result.time_s) + 1) * cell_count - 1);

% Observation summary from the probes, at every accepted step.
summary = readmatrix(files.dam_breach_outputs);
assert(numel(summary) == 33);
maximum_depth = summary(9:13);
assert(max(abs(maximum_depth(1:3) - max(result.probes.depth_m, [], 1))) < 1e-9);
assert(all(isnan(maximum_depth(4:5))));
assert(summary(2) == input.duration_s / 60);

% Without matching probes, the observation points read the stored history.
history_only = hydrohp_postprocess_varying(result, input, ...
    'output_folder', folder, 'run_id', "H", 'figures', false, ...
    'videos', false, 'observation_x_m', [100; 500]);
summary = readmatrix(history_only.dam_breach_outputs);
[~, first_cell] = min(abs(result.x_m - 100));
assert(abs(summary(9) - max(depth(:, first_cell))) < 1e-9);

try
    hydrohp_postprocess_varying(result, input, 'output_folder', folder, ...
        'observation_x_m', (0:100:500).', 'figures', false, 'videos', false);
    error('test:ExpectedFailure', 'Six observation points should fail.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:TooManyObservationPoints'));
end

metrics = struct('file_count', numel(expected), ...
    'stored_times', numel(result.time_s), 'cells', cell_count);
fprintf('HydroHP varying post-processing checks passed: %d files, %d times x %d cells.\n', ...
    metrics.file_count, metrics.stored_times, metrics.cells);
end

function sections = channel(invert_m)
chainage_m = [0; 300; 600];
rows = table();
for index = 1:numel(chainage_m)
    section = table( ...
        repmat("S" + index, 4, 1), repmat(chainage_m(index), 4, 1), ...
        [0; 4; 12; 16], invert_m(index) + [4; 0; 0; 4], ...
        0.035 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section]; %#ok<AGROW>
end
sections = hydrohp_build_section_tables(rows, 0.02, "scm");
end

function remove_folder(folder)
if isfolder(folder)
    rmdir(folder, 's');
end
end
