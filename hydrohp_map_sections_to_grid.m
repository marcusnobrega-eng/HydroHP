function [local_sections, mapping] = hydrohp_map_sections_to_grid(section_data, grid_x_m)
%HYDROHP_MAP_SECTIONS_TO_GRID Map section-property curves to grid locations.

required_fields = {'section_id', 'chainage_m', 'invert_elevation_m', ...
    'maximum_water_surface_elevation_m', 'method', 'table'};
if ~isstruct(section_data) || isempty(section_data) || ...
        ~all(isfield(section_data, required_fields))
    error('hydrohp:InvalidSectionData', ...
        'section_data must come from hydrohp_build_section_tables.');
end
if ~isnumeric(grid_x_m) || ~isvector(grid_x_m) || isempty(grid_x_m) || ...
        any(~isfinite(grid_x_m)) || any(diff(grid_x_m(:)) <= 0)
    error('hydrohp:InvalidGrid', 'grid_x_m must be a finite, strictly increasing vector.');
end

chainage = [section_data.chainage_m].';
if any(~isfinite(chainage)) || any(diff(chainage) <= 0)
    error('hydrohp:InvalidSectionOrder', ...
        'Section chainages must be finite, unique, and strictly increasing.');
end
methods = string({section_data.method});
if numel(unique(methods)) ~= 1
    error('hydrohp:MixedSectionMethods', 'All mapped sections must use the same hydraulic method.');
end

grid_x_m = grid_x_m(:);
tolerance = 1e-10 * max(1, chainage(end) - chainage(1));
if grid_x_m(1) < chainage(1) - tolerance || grid_x_m(end) > chainage(end) + tolerance
    error('hydrohp:GridOutsideSections', ...
        'Grid locations must remain between the first and last section chainages.');
end
grid_x_m(abs(grid_x_m - chainage(1)) <= tolerance) = chainage(1);
grid_x_m(abs(grid_x_m - chainage(end)) <= tolerance) = chainage(end);

count = numel(grid_x_m);
upstream = zeros(count, 1);
downstream = zeros(count, 1);
weight = zeros(count, 1);
template = struct('x_m', 0, 'upstream_section_id', "", ...
    'downstream_section_id', "", 'downstream_weight', 0, 'method', methods(1), ...
    'invert_elevation_m', 0, 'maximum_water_surface_elevation_m', 0, ...
    'table', table());
local_sections = repmat(template, count, 1);

for index = 1:count
    exact = find(abs(chainage - grid_x_m(index)) <= tolerance, 1);
    if ~isempty(exact)
        upstream(index) = exact;
        downstream(index) = exact;
        local_sections(index) = copy_exact(section_data(exact), grid_x_m(index));
        continue
    end

    downstream(index) = find(chainage > grid_x_m(index), 1, 'first');
    upstream(index) = downstream(index) - 1;
    weight(index) = (grid_x_m(index) - chainage(upstream(index))) / ...
        (chainage(downstream(index)) - chainage(upstream(index)));
    local_sections(index) = interpolate_pair(section_data(upstream(index)), ...
        section_data(downstream(index)), grid_x_m(index), weight(index));
end

mapping = struct('grid_x_m', grid_x_m, 'upstream_index', upstream, ...
    'downstream_index', downstream, 'downstream_weight', weight, ...
    'exact_section', upstream == downstream);
end

function local = copy_exact(section, x_m)
local = struct('x_m', x_m, 'upstream_section_id', section.section_id, ...
    'downstream_section_id', section.section_id, 'downstream_weight', 0, ...
    'method', section.method, 'invert_elevation_m', section.invert_elevation_m, ...
    'maximum_water_surface_elevation_m', ...
    section.maximum_water_surface_elevation_m, 'table', section.table);
end

function local = interpolate_pair(upstream, downstream, x_m, weight)
upstream_depth = upstream.maximum_water_surface_elevation_m - upstream.invert_elevation_m;
downstream_depth = downstream.maximum_water_surface_elevation_m - downstream.invert_elevation_m;
relative_depth = uniquetol([upstream.table.depth_m / upstream_depth; ...
    downstream.table.depth_m / downstream_depth], 1e-10, 'DataScale', 1);
relative_depth(1) = 0;
relative_depth(end) = 1;

fields = {'area_m2', 'wetted_perimeter_m', 'top_width_m', ...
    'hydrostatic_integral_m3', 'conveyance_m3s'};
values = zeros(numel(relative_depth), numel(fields));
for index = 1:numel(fields)
    upstream_values = interp1(upstream.table.depth_m, upstream.table.(fields{index}), ...
        relative_depth * upstream_depth, 'pchip');
    downstream_values = interp1(downstream.table.depth_m, downstream.table.(fields{index}), ...
        relative_depth * downstream_depth, 'pchip');
    values(:, index) = (1 - weight) * upstream_values + weight * downstream_values;
end

invert = (1 - weight) * upstream.invert_elevation_m + ...
    weight * downstream.invert_elevation_m;
maximum_wse = (1 - weight) * upstream.maximum_water_surface_elevation_m + ...
    weight * downstream.maximum_water_surface_elevation_m;
depth = relative_depth * (maximum_wse - invert);
water_surface = invert + depth;
area = values(:, 1);
perimeter = values(:, 2);
width = values(:, 3);
hydrostatic = values(:, 4);
conveyance = values(:, 5);
radius = zeros(size(area));
centroid_depth = zeros(size(area));
effective_n = zeros(size(area));
wet = area > 0 & perimeter > 0 & conveyance > 0;
radius(wet) = area(wet) ./ perimeter(wet);
centroid_depth(wet) = hydrostatic(wet) ./ area(wet);
effective_n(wet) = area(wet) .* radius(wet).^(2/3) ./ conveyance(wet);
effective_n(~wet) = effective_n(find(wet, 1, 'first'));

if any(diff(area) <= 0) || any(diff(perimeter) < -1e-12) || ...
        any(~isfinite([depth; water_surface; values(:); radius; centroid_depth; effective_n]))
    error('hydrohp:InvalidMappedSection', ...
        'Interpolated section at x = %.6g m is not finite and monotonic.', x_m);
end

local = struct('x_m', x_m, 'upstream_section_id', upstream.section_id, ...
    'downstream_section_id', downstream.section_id, 'downstream_weight', weight, ...
    'method', upstream.method, 'invert_elevation_m', invert, ...
    'maximum_water_surface_elevation_m', maximum_wse, ...
    'table', table(depth, water_surface, area, perimeter, width, radius, ...
    hydrostatic, centroid_depth, effective_n, conveyance, ...
    'VariableNames', {'depth_m', 'water_surface_elevation_m', 'area_m2', ...
    'wetted_perimeter_m', 'top_width_m', 'hydraulic_radius_m', ...
    'hydrostatic_integral_m3', 'centroid_depth_m', ...
    'effective_manning_n', 'conveyance_m3s'}));
end
