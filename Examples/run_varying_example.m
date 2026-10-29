% Run the varying-geometry kernel on the three example sections, from dry.
repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);

section_data = hydrohp_build_section_tables(fullfile(repository_root, ...
    'Examples', 'cross_sections_example.csv'), 0.02, "scm");

input = struct( ...
    'section_data', section_data, ...
    'face_x_m', linspace(0, 1000, 41).', ...   % 40 cells of 25 m
    'hydrograph_time_s', [0; 600; 1800; 3600; 7200], ...
    'hydrograph_discharge_m3s', [0; 40; 15; 0; 0], ...
    'duration_s', 7200, ...
    'downstream', struct('type', "normal_depth", 'slope', 1e-3), ...
    'initialization_type', "dry", ...
    'output_interval_s', 60, ...
    'probe_x_m', [0; 500; 1000]);
result = hydrohp_run_varying(input);

fprintf('Status: %s\n', result.status);
fprintf('Inflow %.1f m3, outflow %.1f m3, storage %.1f m3, relative residual %.2g\n', ...
    result.ledger.cumulative_inflow_m3, result.ledger.cumulative_outflow_m3, ...
    result.ledger.final_storage_m3, result.ledger.relative_residual);

% Bed and section top at the cell centres, from the tables the solver used.
cell_sections = hydrohp_map_sections_to_grid(section_data, result.x_m);
bed_m = [cell_sections.invert_elevation_m].';
top_m = [cell_sections.maximum_water_surface_elevation_m].';

figure('Name', 'Maximum water level');
plot(result.x_m, top_m, 'k--', result.x_m, result.maxima.stage_m, 'b-', ...
    result.x_m, bed_m, 'k-', 'LineWidth', 1.2);
xlabel('Chainage (m)');
ylabel('Elevation (m)');
legend('Section top', 'Maximum water level', 'Bed', ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
grid on

figure('Name', 'Probe series');
tiledlayout(2, 1);
nexttile
plot(result.probes.time_s / 60, result.probes.face_mass_flux_m3s, 'LineWidth', 1.2);
ylabel('Discharge (m^3/s)');
legend(compose('x = %g m', result.probes.x_m), 'Location', 'northeast');
grid on
nexttile
plot(result.probes.time_s / 60, result.probes.depth_m, 'LineWidth', 1.2);
xlabel('Time (min)');
ylabel('Depth (m)');
grid on
