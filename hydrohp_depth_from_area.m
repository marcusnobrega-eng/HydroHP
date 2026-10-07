function depth_m = hydrohp_depth_from_area(section, area_m2)
%HYDROHP_DEPTH_FROM_AREA Recover depth from one section property table.

if ~isstruct(section) || ~isscalar(section) || ~isfield(section, 'table') || ...
        ~all(ismember({'area_m2', 'depth_m'}, section.table.Properties.VariableNames))
    error('hydrohp:InvalidSectionTable', 'section must contain an H4 hydraulic-property table.');
end
if ~isnumeric(area_m2) || any(~isfinite(area_m2(:))) || any(area_m2(:) < 0)
    error('hydrohp:InvalidAreaQuery', 'area_m2 must contain finite nonnegative values.');
end
maximum_area = section.table.area_m2(end);
if any(area_m2(:) > maximum_area)
    error('hydrohp:AreaOutsideSectionRange', ...
        'Requested area exceeds the valid range for section %s.', section.section_id);
end

depth_m = interp1(section.table.area_m2, section.table.depth_m, ...
    area_m2, 'pchip');
end
