function report = hydrohp_validate_cross_sections(sections)
%HYDROHP_VALIDATE_CROSS_SECTIONS Validate the long-table section schema.

required = {'section_id', 'chainage_m', 'offset_m', 'elevation_m', 'manning_n'};
report = struct('valid', false, 'errors', strings(0, 1), ...
    'warnings', strings(0, 1), 'section_count', 0, 'point_count', 0, ...
    'chainage_m', zeros(0, 1));

if ~istable(sections)
    report.errors(end + 1) = "Cross-sections must be supplied as a MATLAB table.";
    return
end

missing = setdiff(required, sections.Properties.VariableNames, 'stable');
if ~isempty(missing)
    report.errors(end + 1) = "Missing required columns: " + strjoin(missing, ', ') + ".";
    return
end
if isempty(sections)
    report.errors(end + 1) = "The cross-section table is empty.";
    return
end

ids = string(sections.section_id);
if any(ismissing(ids) | strlength(strtrim(ids)) == 0)
    report.errors(end + 1) = "section_id contains missing or empty values.";
end

numeric_columns = {'chainage_m', 'offset_m', 'elevation_m', 'manning_n'};
optional_coordinates = {'x_m', 'y_m'};
has_x = ismember('x_m', sections.Properties.VariableNames);
has_y = ismember('y_m', sections.Properties.VariableNames);
if xor(has_x, has_y)
    report.errors(end + 1) = "Optional coordinates x_m and y_m must be supplied together.";
elseif has_x
    numeric_columns = [numeric_columns, optional_coordinates];
end

for index = 1:numel(numeric_columns)
    name = numeric_columns{index};
    values = sections.(name);
    if ~isnumeric(values) || ~isvector(values) || numel(values) ~= height(sections)
        report.errors(end + 1) = name + " must be one numeric value per row.";
    elseif any(~isfinite(values))
        report.errors(end + 1) = name + " contains nonfinite values.";
    end
end
if ~isempty(report.errors)
    return
end

if any(sections.chainage_m < 0)
    report.errors(end + 1) = "chainage_m must be nonnegative.";
end
if any(sections.manning_n <= 0)
    report.errors(end + 1) = "manning_n must be positive.";
end

section_ids = unique(ids, 'stable');
chainages = nan(numel(section_ids), 1);
for index = 1:numel(section_ids)
    rows = ids == section_ids(index);
    section_chainage = unique(sections.chainage_m(rows));
    offsets = sections.offset_m(rows);
    if numel(section_chainage) ~= 1
        report.errors(end + 1) = "Section " + section_ids(index) + ...
            " must have one chainage_m value.";
    else
        chainages(index) = section_chainage;
    end
    if nnz(rows) < 3
        report.errors(end + 1) = "Section " + section_ids(index) + ...
            " must contain at least three profile points.";
    end
    if numel(unique(offsets)) ~= numel(offsets)
        report.errors(end + 1) = "Section " + section_ids(index) + ...
            " contains duplicate offset_m values.";
    end
end

finite_chainages = chainages(isfinite(chainages));
if numel(unique(finite_chainages)) ~= numel(finite_chainages)
    report.errors(end + 1) = "Each section_id must have a unique chainage_m.";
end

ordered = sortrows(table(sections.chainage_m, ids, sections.offset_m), [1, 2, 3]);
current = table(sections.chainage_m, ids, sections.offset_m);
if ~isequal(current, ordered)
    report.warnings(end + 1) = ...
        "Rows are not ordered by chainage_m, section_id, and offset_m; the reader will reorder them.";
end

report.section_count = numel(section_ids);
report.point_count = height(sections);
report.chainage_m = sort(finite_chainages);
report.valid = isempty(report.errors);
end
