function metrics = test_varying_run
%TEST_VARYING_RUN Check H6 boundaries, initialization, and mass ledger.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
sections = sloping_sections;
base = struct( ...
    'section_data', sections, ...
    'face_x_m', linspace(0, 100, 31).', ...
    'duration_s', 20, ...
    'downstream', struct('type', "normal_depth", 'slope', 0.001), ...
    'base_discharge_m3s', 1, ...
    'maximum_dt_s', 0.2, ...
    'spinup_duration_s', 100, ...
    'steady_tolerance_per_s', 1e-4);

% Constant base flow exercises normal-depth initialization and both boundaries.
input = base;
input.hydrograph_time_s = [0; 20];
input.hydrograph_discharge_m3s = [1; 1];
constant = hydrohp_run_varying(input);
assert(constant.initialization.converged);
assert(constant.initialization.relative_state_change_rate_per_s < 1e-4);
assert(constant.initialization.relative_discharge_error < 0.03);
assert(abs(constant.ledger.relative_residual) < 1e-11);
assert(all(isfinite([constant.area_m2(:); constant.discharge_m3s(:)])));

% A piecewise-linear hydrograph must enter with its exact integrated volume.
input = base;
input.spinup_duration_s = 0;
input.hydrograph_time_s = [0; 10; 20];
input.hydrograph_discharge_m3s = [1; 1.5; 1];
event = hydrohp_run_varying(input);
assert(abs(event.ledger.cumulative_inflow_m3 - 25) < 1e-10);
assert(abs(event.ledger.relative_residual) < 1e-11);

% A user stage boundary must preserve a motionless horizontal surface.
input = base;
input.downstream = struct('type', "stage", 'stage_m', 0.8);
input.initialization_type = "stage";
input.initial_stage_m = 0.8;
input.base_discharge_m3s = 0;
input.spinup_duration_s = 0;
input.include_friction = false;
input.hydrograph_time_s = [0; 20];
input.hydrograph_discharge_m3s = [0; 0];
still = hydrohp_run_varying(input);
assert(max(abs(still.discharge_m3s(:))) < 1e-9);
assert(abs(still.ledger.relative_residual) < 1e-11);

metrics = struct( ...
    'steady_relative_discharge_error', ...
        constant.initialization.relative_discharge_error, ...
    'steady_relative_state_change_rate_per_s', ...
        constant.initialization.relative_state_change_rate_per_s, ...
    'steady_relative_mass_residual', constant.ledger.relative_residual, ...
    'hydrograph_cumulative_inflow_m3', event.ledger.cumulative_inflow_m3, ...
    'hydrograph_relative_mass_residual', event.ledger.relative_residual, ...
    'still_maximum_absolute_discharge_m3s', ...
        max(abs(still.discharge_m3s(:))), ...
    'still_relative_mass_residual', still.ledger.relative_residual);
fprintf('HydroHP varying-section run checks passed.\n');
end

function sections = sloping_sections
chainage_m = [0; 50; 100];
invert_m = [0.1; 0.05; 0];
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
