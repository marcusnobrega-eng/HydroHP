function metrics = test_overtopping_guard
%TEST_OVERTOPPING_GUARD A flow beyond the section tables is an applicability
% failure. With overtopping_action = "stop" the run ends at the last accepted
% state and reports where and when; the default still raises the error.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
input = struct( ...
    'section_data', shallow_sections(0.8), ...
    'face_x_m', linspace(0, 300, 31).', ...
    'hydrograph_time_s', [0; 60; 600], ...
    'hydrograph_discharge_m3s', [0.5; 40; 40], ...
    'duration_s', 600, ...
    'downstream', struct('type', "normal_depth", 'slope', 0.002), ...
    'base_discharge_m3s', 0.5, ...
    'spinup_duration_s', 120, ...
    'steady_tolerance_per_s', 1e9, ...
    'maximum_dt_s', 1, ...
    'reconstruction', "hydrostatic", ...
    'output_interval_s', 10);

identifier = "";
try
    hydrohp_run_varying(input);
catch exception
    identifier = string(exception.identifier);
end
assert(ismember(identifier, ["hydrohp:AreaOutsideSectionRange", ...
    "hydrohp:StageOutsideSectionRange"]), ...
    'Expected an out-of-range error by default, got "%s".', identifier);

input.overtopping_action = "stop";
result = hydrohp_run_varying(input);
event = result.overtopping;
assert(result.status == "overtopped");
assert(event.phase == "event" && event.identifier == identifier);
assert(event.model_time_s > 0 && event.model_time_s < input.duration_s);
assert(isfinite(event.x_m) && event.x_m >= 0 && event.x_m <= 300);
assert(result.time_s(end) == event.model_time_s);
assert(result.ledger.step_count == event.accepted_step_count);
assert(abs(result.ledger.relative_residual) < 1e-12);
assert(all(isfinite(result.maxima.depth_m)));
assert(result.metadata.overtopping_action == "stop");

% A run that stays inside its tables reports completion.
input.section_data = shallow_sections(4);
contained = hydrohp_run_varying(input);
assert(contained.status == "completed" && isempty(contained.overtopping));

metrics = struct('overtopping_time_s', event.model_time_s, ...
    'overtopping_x_m', event.x_m);
fprintf('HydroHP overtopping-guard checks passed: stop at t = %.1f s, x = %.1f m.\n', ...
    event.model_time_s, event.x_m);
end

function sections = shallow_sections(bank_height_m)
chainage_m = [0; 150; 300];
rows = table();
for index = 1:numel(chainage_m)
    invert = 1 - 0.002 * chainage_m(index);
    section = table( ...
        repmat("S" + index, 4, 1), repmat(chainage_m(index), 4, 1), ...
        [0; 1; 9; 10], invert + bank_height_m * [1; 0; 0; 1], ...
        0.035 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section]; %#ok<AGROW>
end
sections = hydrohp_build_section_tables(rows, 0.02, "scm");
end
