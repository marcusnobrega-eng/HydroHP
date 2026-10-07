function input = hydrohp_read_excel(excel_file)
%HYDROHP_READ_EXCEL Create a HydroHP input object from a legacy workbook.

if ~(ischar(excel_file) || (isstring(excel_file) && isscalar(excel_file)))
    error('hydrohp:InvalidInputFile', 'excel_file must be a character vector or scalar string.');
end

[exists, attributes] = fileattrib(char(excel_file));
if ~exists || attributes.directory
    error('hydrohp:InputFileNotFound', 'HydroHP input workbook not found: %s', char(excel_file));
end

input = struct();
input.excel_file = attributes.Name;
input.postprocess = false;
input.verbose = false;
end
