function metrics = test_pool_initialization
%TEST_POOL_INITIALIZATION A bed with a pool (locally rising invert) starts
% from normal depth at the reach-average slope in the non-positive cells, and
% the spin-up still reaches a steady state without creating volume. A reach
% that rises overall cannot be initialized at normal depth.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
input = struct( ...
    'section_data', sections([0.10; 0.14; 0]), ...
    'face_x_m', linspace(0, 100, 31).', ...
    'hydrograph_time_s', [0; 20], ...
    'hydrograph_discharge_m3s', [1; 1], ...
    'duration_s', 20, ...
    'downstream', struct('type', "normal_depth", 'slope', 0.001), ...
    'base_discharge_m3s', 1, ...
    'maximum_dt_s', 0.2, ...
    'spinup_duration_s', 300, ...
    'steady_tolerance_per_s', 1e-4);
result = hydrohp_run_varying(input);
slope = result.initialization.initial_slope;
assert(slope.nonpositive_cell_count > 0 && slope.nonpositive_cell_count < 30);
assert(abs(slope.reach_average_slope - 0.001) < 1e-12);
assert(result.initialization.converged);
assert(result.initialization.ledger.wet_dry_adjustment_m3 == 0);
assert(result.status == "completed");
assert(abs(result.ledger.relative_residual) < 1e-11);

% A downhill bed records no fallback cells.
input.section_data = sections([0.10; 0.05; 0]);
downhill = hydrohp_run_varying(input);
assert(downhill.initialization.initial_slope.nonpositive_cell_count == 0);

% A reach that rises overall has no normal-depth starting state.
input.section_data = sections([0; 0.05; 0.10]);
input.downstream.slope = 0.001;
try
    hydrohp_run_varying(input);
    error('hydrohp:TestFailed', 'An adverse reach was initialized at normal depth.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:InvalidInitialization'));
end

metrics = struct('pool_fallback_cells', slope.nonpositive_cell_count, ...
    'spinup_change_rate_per_s', result.initialization.relative_state_change_rate_per_s);
fprintf('HydroHP pool-initialization checks passed: %d fallback cells, converged.\n', ...
    slope.nonpositive_cell_count);
end

function data = sections(invert_m)
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
data = hydrohp_build_section_tables(rows, 0.01, "scm");
end
