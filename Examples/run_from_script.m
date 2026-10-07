% Run HydroHP without editing the solver script.
repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);

input = hydrohp_read_excel(fullfile(repository_root, 'Examples', ...
    'Example_2_Dam_Break.xlsx'));

% Optional script-defined values:
% input.hydrograph.time_min = [0; 5; 15; 30];
% input.hydrograph.discharge_m3s = [5; 100; 40; 5];
% input.manning_n = 0.035;

result = hydrohp_run(input);
disp(result)
