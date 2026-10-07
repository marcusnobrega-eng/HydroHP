function test_cross_section_input
%TEST_CROSS_SECTION_INPUT Check table, CSV, Excel, ordering, and QA failures.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
csv_file = fullfile(repository_root, 'Examples', 'cross_sections_example.csv');
xlsx_file = fullfile(repository_root, 'Examples', 'HydroHP_MultiSection_Input.xlsx');

[from_csv, csv_report] = hydrohp_read_cross_sections(csv_file);
[from_xlsx, xlsx_report] = hydrohp_read_cross_sections(xlsx_file);
assert(csv_report.valid && xlsx_report.valid);
assert(csv_report.section_count == 3 && csv_report.point_count == 15);
assert(isequal(from_csv, from_xlsx));

shuffled = from_csv([6:10, 1:5, 11:15], :);
[from_table, table_report] = hydrohp_read_cross_sections(shuffled);
assert(table_report.valid && ~isempty(table_report.warnings));
assert(isequal(from_csv, from_table));

invalid = from_csv;
invalid.chainage_m(invalid.section_id == "XS002") = 0;
report = hydrohp_validate_cross_sections(invalid);
assert(~report.valid && any(contains(report.errors, "unique chainage_m")));

invalid = from_csv;
invalid.offset_m(2) = invalid.offset_m(1);
report = hydrohp_validate_cross_sections(invalid);
assert(~report.valid && any(contains(report.errors, "duplicate offset_m")));

invalid = removevars(from_csv, 'manning_n');
report = hydrohp_validate_cross_sections(invalid);
assert(~report.valid && any(contains(report.errors, "Missing required columns")));

invalid = from_csv;
invalid.manning_n(1) = 0;
report = hydrohp_validate_cross_sections(invalid);
assert(~report.valid && any(contains(report.errors, "must be positive")));

invalid = from_csv(from_csv.section_id ~= "XS003" | from_csv.offset_m < 20, :);
report = hydrohp_validate_cross_sections(invalid);
assert(~report.valid && any(contains(report.errors, "at least three")));

fprintf('HydroHP cross-section schema checks passed.\n');
end
