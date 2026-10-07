function result = hydrohp_run(input)
%HYDROHP_RUN Run the legacy HydroHP solver through a callable interface.
%
% result = hydrohp_run("Examples/Example_2_Dam_Break.xlsx")
% result = hydrohp_run(hydrohp_read_excel("Examples/Example_2_Dam_Break.xlsx"))
%
% A script-created input may optionally define:
%   input.hydrograph.time_min
%   input.hydrograph.discharge_m3s
%   input.manning_n

if ischar(input) || (isstring(input) && isscalar(input))
    input = hydrohp_read_excel(input);
elseif ~isstruct(input) || ~isscalar(input)
    error('hydrohp:InvalidInput', 'input must be an Excel path or scalar structure.');
end

if ~isfield(input, 'excel_file')
    error('hydrohp:MissingInputFile', 'input.excel_file is required.');
end
input = apply_defaults(input);
validate_input(input);

[exists, attributes] = fileattrib(char(input.excel_file));
if ~exists || attributes.directory
    error('hydrohp:InputFileNotFound', 'HydroHP input workbook not found: %s', char(input.excel_file));
end

HydroHP_Input_File = attributes.Name;
HydroHP_Input_Overrides = input;
HydroHP_Options = struct( ...
    'postprocess', input.postprocess, ...
    'clear_command_window', false);

warning_state = warning;
HydroHP_Cleanup = onCleanup(@() warning(warning_state));

if input.verbose
    HydroHP_1D_V5;
    console_output = '';
else
    console_output = evalc('HydroHP_1D_V5');
end

result = struct();
result.input_file = HydroHP_Input_File;
result.time_s = time_save(:);
result.x_m = x(:);
result.flow_area_m2 = Flow_Area;
result.discharge_m3s = Discharge;
result.depth_m = Depth;
result.velocity_ms = Velocity;
result.froude = Froude;
result.courant = Courant;
result.execution_time_s = computational_time;
result.console_output = console_output;
end

function input = apply_defaults(input)
if ~isfield(input, 'postprocess')
    input.postprocess = false;
end
if ~isfield(input, 'verbose')
    input.verbose = false;
end
end

function validate_input(input)
if ~(ischar(input.excel_file) || (isstring(input.excel_file) && isscalar(input.excel_file)))
    error('hydrohp:InvalidInputFile', 'input.excel_file must be a character vector or scalar string.');
end
if ~(islogical(input.postprocess) && isscalar(input.postprocess))
    error('hydrohp:InvalidPostprocess', 'input.postprocess must be a logical scalar.');
end
if ~(islogical(input.verbose) && isscalar(input.verbose))
    error('hydrohp:InvalidVerbose', 'input.verbose must be a logical scalar.');
end

if isfield(input, 'hydrograph') && ~isempty(input.hydrograph)
    hydrograph = input.hydrograph;
    if ~isstruct(hydrograph) || ~isscalar(hydrograph) || ...
            ~isfield(hydrograph, 'time_min') || ~isfield(hydrograph, 'discharge_m3s')
        error('hydrohp:InvalidHydrograph', ...
            'input.hydrograph requires time_min and discharge_m3s vectors.');
    end
    time_min = hydrograph.time_min;
    discharge = hydrograph.discharge_m3s;
    valid = isnumeric(time_min) && isvector(time_min) && numel(time_min) >= 2 && ...
        isnumeric(discharge) && isvector(discharge) && numel(discharge) == numel(time_min) && ...
        all(isfinite(time_min)) && all(isfinite(discharge)) && ...
        time_min(1) == 0 && all(diff(time_min) > 0) && all(discharge >= 0);
    if ~valid
        error('hydrohp:InvalidHydrograph', ...
            'Hydrograph vectors must be finite, nonnegative, equal-length, start at zero, and have increasing time.');
    end
end

if isfield(input, 'manning_n') && ~isempty(input.manning_n)
    if ~(isnumeric(input.manning_n) && isscalar(input.manning_n) && ...
            isfinite(input.manning_n) && input.manning_n > 0)
        error('hydrohp:InvalidManning', 'input.manning_n must be a positive finite scalar.');
    end
end
end
