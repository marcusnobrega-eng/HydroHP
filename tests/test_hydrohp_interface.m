function test_hydrohp_interface
%TEST_HYDROHP_INTERFACE Check equivalent Excel-path and structure calls.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
example_file = fullfile(repository_root, 'Examples', 'Example_2_Dam_Break.xlsx');

% Read the legacy hydrograph once, then supply the same values explicitly as
% script data to exercise the programmatic override path.
HydroHP_Input_File = example_file;
Read_Input_Data;
script_hydrograph = struct('time_min', time, 'discharge_m3s', Qe1);

from_path = hydrohp_run(example_file);
script_input = hydrohp_read_excel(example_file);
script_input.hydrograph = script_hydrograph;
from_struct = hydrohp_run(script_input);

fields = {'time_s', 'x_m', 'flow_area_m2', 'discharge_m3s', ...
    'depth_m', 'velocity_ms', 'froude', 'courant'};
for index = 1:numel(fields)
    field = fields{index};
    assert(isequaln(from_path.(field), from_struct.(field)), ...
        'Excel-path and structure results differ for %s.', field);
end

invalid = hydrohp_read_excel(example_file);
invalid.hydrograph.time_min = [0; 2; 1];
invalid.hydrograph.discharge_m3s = [1; 2; 3];
try
    hydrohp_run(invalid);
    error('hydrohp:TestFailed', 'Invalid hydrograph was accepted.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:InvalidHydrograph'));
end

fprintf('HydroHP callable-interface checks passed.\n');
end
