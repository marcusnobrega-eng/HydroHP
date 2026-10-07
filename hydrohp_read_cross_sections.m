function [sections, report] = hydrohp_read_cross_sections(source)
%HYDROHP_READ_CROSS_SECTIONS Read and validate multi-section profile data.
%
% source may be a MATLAB table, a CSV file, or an Excel workbook containing
% a sheet named CrossSections.

if istable(source)
    sections = source;
elseif ischar(source) || (isstring(source) && isscalar(source))
    source = char(source);
    if ~isfile(source)
        error('hydrohp:CrossSectionFileNotFound', ...
            'Cross-section input not found: %s', source);
    end
    [~, ~, extension] = fileparts(source);
    switch lower(extension)
        case '.csv'
            sections = readtable(source, 'VariableNamingRule', 'preserve');
        case {'.xlsx', '.xls'}
            sections = readtable(source, 'Sheet', 'CrossSections', ...
                'VariableNamingRule', 'preserve');
        otherwise
            error('hydrohp:UnsupportedCrossSectionFile', ...
                'Cross-sections must use CSV, XLSX, or XLS format.');
    end
else
    error('hydrohp:InvalidCrossSectionSource', ...
        'source must be a MATLAB table or a CSV/Excel path.');
end

report = hydrohp_validate_cross_sections(sections);
if ~report.valid
    error('hydrohp:InvalidCrossSections', '%s', ...
        strjoin(cellstr(report.errors), newline));
end

sections.section_id = string(sections.section_id);
sections = sortrows(sections, {'chainage_m', 'section_id', 'offset_m'});
end
