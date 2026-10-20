function metrics = test_event_outputs
%TEST_EVENT_OUTPUTS P3: sparse history, online maxima, and probe series.
%
% Sparse output must not change the solution: the final state, ledger,
% maxima, and probes equal those of a run that stores every step, and the
% sparse history rows are rows of the full history. Online maxima must equal
% maxima recomputed from the full history.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
input = struct( ...
    'section_data', trapezoid_sections([0.1; 0.05; 0]), ...
    'face_x_m', linspace(0, 100, 31).', ...
    'hydrograph_time_s', [0; 5; 10; 30], ...
    'hydrograph_discharge_m3s', [1; 4; 1.5; 1], ...
    'duration_s', 30, ...
    'downstream', struct('type', "normal_depth", 'slope', 0.001), ...
    'base_discharge_m3s', 1, ...
    'maximum_dt_s', 0.02, ...
    'spinup_duration_s', 50, ...
    'steady_tolerance_per_s', 1e9, ...
    'reconstruction', "hydrostatic", ...
    'probe_x_m', [10; 55; 99]);
full = hydrohp_run_varying(input);
input.output_interval_s = 2.5;
sparse = hydrohp_run_varying(input);

assert(isequal(sparse.final_area_m2, full.final_area_m2));
assert(isequal(sparse.final_discharge_m3s, full.final_discharge_m3s));
assert(isequaln(sparse.ledger, full.ledger));
assert(isequaln(sparse.maxima, full.maxima));
assert(isequaln(sparse.probes, full.probes));
assert(numel(full.time_s) == full.ledger.step_count + 1);
assert(numel(full.time_s) > 1024, 'The test must exercise buffer growth.');
assert(numel(sparse.time_s) < numel(full.time_s) / 5);
[present, rows] = ismember(sparse.time_s, full.time_s);
assert(all(present));
assert(isequal(sparse.area_m2, full.area_m2(rows, :)));
assert(isequal(sparse.discharge_m3s, full.discharge_m3s(rows, :)));
assert(sparse.time_s(end) == input.duration_s);
assert(all(diff(floor(sparse.time_s(2:end-1) / 2.5)) >= 1));

% Online maxima equal maxima over every stored step.
cells = hydrohp_map_sections_to_grid(input.section_data, full.x_m);
depth = zeros(size(full.area_m2));
for index = 1:numel(cells)
    depth(:, index) = hydrohp_depth_from_area(cells(index), full.area_m2(:, index));
end
speed = abs(full.discharge_m3s) ./ full.area_m2;
speed(full.area_m2 <= 1e-10) = 0;
assert(max(abs(max(depth, [], 1).' - full.maxima.depth_m)) < 1e-12);
assert(isequal(max(full.area_m2, [], 1).', full.maxima.area_m2));
assert(isequal(max(abs(full.discharge_m3s), [], 1).', full.maxima.discharge_abs_m3s));
assert(isequal(max(speed, [], 1).', full.maxima.velocity_m_s));
assert(max(abs(max(depth .* speed, [], 1).' - full.maxima.depth_velocity_m2_s)) < 1e-12);
[~, peak_row] = max(depth, [], 1);
assert(max(abs(full.time_s(peak_row) - full.maxima.time_of_maximum_depth_s)) < 1e-12);

% Probes record every step at the nearest cell and face.
probes = full.probes;
assert(isequal(probes.time_s, full.time_s));
assert(isequal(probes.discharge_m3s, full.discharge_m3s(:, probes.cell_index)));
assert(all(isnan(probes.face_mass_flux_m3s(1, :))));
assert(all(isfinite(probes.face_mass_flux_m3s(2:end, :)), 'all'));
assert(numel(full.maxima.face_mass_flux_abs_m3s) == numel(full.x_m) + 1);

metrics = struct('full_rows', numel(full.time_s), ...
    'sparse_rows', numel(sparse.time_s), ...
    'maximum_depth_m', max(full.maxima.depth_m), ...
    'maximum_velocity_m_s', max(full.maxima.velocity_m_s));
fprintf('HydroHP event-output checks passed: %d full rows, %d sparse rows.\n', ...
    metrics.full_rows, metrics.sparse_rows);
end

function sections = trapezoid_sections(invert_m)
chainage_m = [0; 50; 100];
rows = table();
for index = 1:numel(chainage_m)
    section = table( ...
        repmat("S" + index, 4, 1), repmat(chainage_m(index), 4, 1), ...
        [0; 0.8; 7.2; 8], invert_m(index) + [3; 0; 0; 3], ...
        0.03 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section]; %#ok<AGROW>
end
sections = hydrohp_build_section_tables(rows, 0.01, "scm");
end
