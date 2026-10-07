function test_section_hydraulics
%TEST_SECTION_HYDRAULICS Check analytic geometry and multi-section tables.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);

trapezoid = table( ...
    repmat("TRAP", 4, 1), zeros(4, 1), [0; 2; 8; 10], [2; 0; 0; 2], ...
    0.03 * ones(4, 1), 'VariableNames', {'section_id', 'chainage_m', ...
    'offset_m', 'elevation_m', 'manning_n'});
[scm, report] = hydrohp_build_section_tables(trapezoid, 0.25, "scm");
assert(report.valid && report.section_count == 1);
row = find(scm.table.depth_m == 1, 1);
assert(abs(scm.table.area_m2(row) - 7) < 1e-12);
assert(abs(scm.table.top_width_m(row) - 8) < 1e-12);
assert(abs(scm.table.wetted_perimeter_m(row) - (6 + 2*sqrt(2))) < 1e-12);
assert(abs(scm.table.hydrostatic_integral_m3(row) - (10/3)) < 1e-12);
assert(abs(scm.table.centroid_depth_m(row) - (10/21)) < 1e-12);
assert(max(abs(scm.table.effective_manning_n - 0.03)) < 1e-12);
assert(abs(hydrohp_depth_from_area(scm, 7) - 1) < 1e-12);
assert(all(diff(scm.table.area_m2) > 0));
assert(all(diff(scm.table.wetted_perimeter_m) >= 0));

dcm = hydrohp_build_section_tables(trapezoid, 0.25, "dcm");
assert(max(abs(scm.table.conveyance_m3s - dcm.table.conveyance_m3s)) < 1e-12);

example_file = fullfile(repository_root, 'Examples', 'cross_sections_example.csv');
[sections, multi_report] = hydrohp_build_section_tables(example_file, 0.10, "scm");
assert(multi_report.valid && numel(sections) == 3);
for index = 1:numel(sections)
    values = sections(index).table;
    assert(all(diff(values.area_m2) > 0));
    assert(all(diff(values.wetted_perimeter_m) >= 0));
    assert(all(isfinite(values{:,:}), 'all'));
    assert(values.water_surface_elevation_m(end) == ...
        sections(index).maximum_water_surface_elevation_m);
end

dcm_sections = hydrohp_build_section_tables(example_file, 0.10, "dcm");
for index = 1:numel(dcm_sections)
    assert(all(isfinite(dcm_sections(index).table{:,:}), 'all'));
    assert(all(dcm_sections(index).table.conveyance_m3s(2:end) > 0));
end

try
    hydrohp_depth_from_area(sections(1), sections(1).table.area_m2(end) + 1);
    error('hydrohp:TestFailed', 'Out-of-range area was accepted.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:AreaOutsideSectionRange'));
end

fprintf('HydroHP section-property checks passed.\n');
end
