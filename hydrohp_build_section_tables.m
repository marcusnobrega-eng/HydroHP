function [section_data, report] = hydrohp_build_section_tables(source, depth_step_m, method)
%HYDROHP_BUILD_SECTION_TABLES Build hydraulic-property curves per profile.

if nargin < 2 || isempty(depth_step_m)
    depth_step_m = 0.05;
end
if nargin < 3 || isempty(method)
    method = "scm";
end
if ~(isnumeric(depth_step_m) && isscalar(depth_step_m) && ...
        isfinite(depth_step_m) && depth_step_m > 0)
    error('hydrohp:InvalidDepthStep', 'depth_step_m must be a positive finite scalar.');
end
method = lower(string(method));
if ~isscalar(method) || ~ismember(method, ["scm", "dcm"])
    error('hydrohp:InvalidSectionMethod', 'method must be "scm" or "dcm".');
end

[sections, input_report] = hydrohp_read_cross_sections(source);
section_ids = unique(sections.section_id, 'stable');
template = struct('section_id', "", 'chainage_m', 0, 'points', table(), ...
    'invert_elevation_m', 0, 'maximum_water_surface_elevation_m', 0, ...
    'method', method, 'table', table());
section_data = repmat(template, numel(section_ids), 1);

for index = 1:numel(section_ids)
    points = sections(sections.section_id == section_ids(index), :);
    x = points.offset_m;
    z = points.elevation_m;
    maximum_wse = min(z([1, end]));
    invert = min(z);
    maximum_depth = maximum_wse - invert;
    if maximum_depth <= 0
        error('hydrohp:InvalidSectionGeometry', ...
            'Section %s has no closed stage range above its invert.', section_ids(index));
    end

    vertex_depths = z - invert;
    vertex_depths = vertex_depths(vertex_depths >= 0 & vertex_depths <= maximum_depth);
    depth = uniquetol([0:depth_step_m:maximum_depth, maximum_depth, vertex_depths.'], ...
        1e-10, 'DataScale', 1).';
    depth(end) = maximum_depth;
    count = numel(depth);
    area = zeros(count, 1);
    perimeter = zeros(count, 1);
    width = zeros(count, 1);
    hydrostatic = zeros(count, 1);
    effective_n = zeros(count, 1);
    conveyance = zeros(count, 1);

    for row = 2:count
        water_surface = invert + depth(row);
        [area(row), perimeter(row), width(row), hydrostatic(row), ...
            effective_n(row), conveyance(row)] = section_properties( ...
            x, z, points.manning_n(1:end-1), water_surface, method);
    end

    first_wet = find(effective_n > 0, 1, 'first');
    effective_n(1:first_wet-1) = effective_n(first_wet);
    radius = zeros(count, 1);
    wet = perimeter > 0;
    radius(wet) = area(wet) ./ perimeter(wet);
    centroid_depth = zeros(count, 1);
    centroid_depth(wet) = hydrostatic(wet) ./ area(wet);
    water_surface = invert + depth;

    if any(diff(area) <= 0) || any(diff(perimeter) < -1e-12) || ...
            any(~isfinite([area; perimeter; width; hydrostatic; effective_n; conveyance]))
        error('hydrohp:InvalidHydraulicTable', ...
            'Section %s did not produce finite monotonic area and perimeter curves.', section_ids(index));
    end

    section_data(index).section_id = section_ids(index);
    section_data(index).chainage_m = points.chainage_m(1);
    section_data(index).points = points;
    section_data(index).invert_elevation_m = invert;
    section_data(index).maximum_water_surface_elevation_m = maximum_wse;
    section_data(index).method = method;
    section_data(index).table = table(depth, water_surface, area, perimeter, ...
        width, radius, hydrostatic, centroid_depth, effective_n, conveyance, ...
        'VariableNames', {'depth_m', 'water_surface_elevation_m', 'area_m2', ...
        'wetted_perimeter_m', 'top_width_m', 'hydraulic_radius_m', ...
        'hydrostatic_integral_m3', 'centroid_depth_m', ...
        'effective_manning_n', 'conveyance_m3s'});
end

report = struct('valid', true, 'method', method, ...
    'section_count', numel(section_data), ...
    'input_report', input_report, ...
    'table_rows', arrayfun(@(item) height(item.table), section_data));
end

function [area, perimeter, width, hydrostatic, effective_n, conveyance] = ...
        section_properties(x, z, segment_n, water_surface, method)
segment_count = numel(x) - 1;
segment_area = zeros(segment_count, 1);
segment_perimeter = zeros(segment_count, 1);
segment_width = zeros(segment_count, 1);
segment_hydrostatic = zeros(segment_count, 1);

for index = 1:segment_count
    dx = x(index + 1) - x(index);
    dz = z(index + 1) - z(index);
    depths = water_surface - z(index:index + 1);
    if max(depths) <= 0
        continue
    end
    segment_length = hypot(dx, dz);
    if min(depths) >= 0
        wet_dx = dx;
        wet_length = segment_length;
        d1 = depths(1);
        d2 = depths(2);
    else
        positive_depth = max(depths);
        wet_fraction = positive_depth / abs(diff(depths));
        wet_dx = dx * wet_fraction;
        wet_length = segment_length * wet_fraction;
        d1 = positive_depth;
        d2 = 0;
    end
    segment_area(index) = wet_dx * (d1 + d2) / 2;
    segment_perimeter(index) = wet_length;
    segment_width(index) = wet_dx;
    segment_hydrostatic(index) = wet_dx * (d1^2 + d1*d2 + d2^2) / 6;
end

area = sum(segment_area);
perimeter = sum(segment_perimeter);
width = sum(segment_width);
hydrostatic = sum(segment_hydrostatic);

if method == "scm"
    effective_n = (sum(segment_perimeter .* segment_n.^(3/2)) / perimeter)^(2/3);
    conveyance = area * (area / perimeter)^(2/3) / effective_n;
else
    group_start = [1; find(diff(segment_n) ~= 0) + 1];
    group_end = [group_start(2:end) - 1; segment_count];
    conveyance = 0;
    for index = 1:numel(group_start)
        rows = group_start(index):group_end(index);
        group_area = sum(segment_area(rows));
        group_perimeter = sum(segment_perimeter(rows));
        if group_area > 0 && group_perimeter > 0
            conveyance = conveyance + group_area * ...
                (group_area / group_perimeter)^(2/3) / segment_n(group_start(index));
        end
    end
    effective_n = area * (area / perimeter)^(2/3) / conveyance;
end
end
