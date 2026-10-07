function test_section_grid_mapping
%TEST_SECTION_GRID_MAPPING Check bracketing and property interpolation.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
example_file = fullfile(repository_root, 'Examples', 'cross_sections_example.csv');
sections = hydrohp_build_section_tables(example_file, 0.10, "scm");
grid_x = [0; 250; 500; 750; 1000];
[local, mapping] = hydrohp_map_sections_to_grid(sections, grid_x);

assert(isequal(mapping.upstream_index, [1; 1; 2; 2; 3]));
assert(isequal(mapping.downstream_index, [1; 2; 2; 3; 3]));
assert(isequal(mapping.downstream_weight, [0; 0.5; 0; 0.5; 0]));
assert(isequal(mapping.exact_section, logical([1; 0; 1; 0; 1])));
assert(isequal(local(1).table, sections(1).table));
assert(isequal(local(3).table, sections(2).table));
assert(isequal(local(5).table, sections(3).table));

midpoint = local(2);
relative_depth = midpoint.table.depth_m / ...
    (midpoint.maximum_water_surface_elevation_m - midpoint.invert_elevation_m);
row = find(abs(relative_depth - 0.5) < 1e-12, 1);
upstream_depth = sections(1).maximum_water_surface_elevation_m - ...
    sections(1).invert_elevation_m;
downstream_depth = sections(2).maximum_water_surface_elevation_m - ...
    sections(2).invert_elevation_m;
expected_area = 0.5 * interp1(sections(1).table.depth_m, ...
    sections(1).table.area_m2, 0.5 * upstream_depth, 'pchip') + ...
    0.5 * interp1(sections(2).table.depth_m, ...
    sections(2).table.area_m2, 0.5 * downstream_depth, 'pchip');
assert(abs(midpoint.table.area_m2(row) - expected_area) < 1e-12);
assert(midpoint.invert_elevation_m == 0.5 * ...
    (sections(1).invert_elevation_m + sections(2).invert_elevation_m));
assert(midpoint.maximum_water_surface_elevation_m == 0.5 * ...
    (sections(1).maximum_water_surface_elevation_m + ...
    sections(2).maximum_water_surface_elevation_m));
wet = midpoint.table.area_m2 > 0;
assert(max(abs(midpoint.table.hydraulic_radius_m(wet) - ...
    midpoint.table.area_m2(wet) ./ midpoint.table.wetted_perimeter_m(wet))) < 1e-12);
assert(max(abs(midpoint.table.centroid_depth_m(wet) - ...
    midpoint.table.hydrostatic_integral_m3(wet) ./ midpoint.table.area_m2(wet))) < 1e-12);

try
    hydrohp_map_sections_to_grid(sections, [-1; 0; 500]);
    error('hydrohp:TestFailed', 'Out-of-range grid was accepted.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:GridOutsideSections'));
end

fprintf('HydroHP section-to-grid mapping checks passed.\n');
end
