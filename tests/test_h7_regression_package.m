function metrics = test_h7_regression_package
%TEST_H7_REGRESSION_PACKAGE Published normal flow and analytical dam break.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);

% Published triangular-section values from Gomes Jr. et al. (2023), Table 1.
triangle = table(repmat("T", 3, 1), zeros(3, 1), [0; 1; 3], [2; 0; 2], ...
    0.02 * ones(3, 1), 'VariableNames', {'section_id', 'chainage_m', ...
    'offset_m', 'elevation_m', 'manning_n'});
section = hydrohp_build_section_tables(triangle, 0.01, "scm");
row = height(section.table);
area_m2 = section.table.area_m2(row);
perimeter_m = section.table.wetted_perimeter_m(row);
radius_m = section.table.hydraulic_radius_m(row);
conveyance_m3s = section.table.conveyance_m3s(row);
discharge_m3s = conveyance_m3s * sqrt(0.001);
velocity_m_s = discharge_m3s / area_m2;
published = [5.06, 0.59, 105.58, 3.339, 1.11];
computed = [perimeter_m, radius_m, conveyance_m3s, discharge_m3s, velocity_m_s];
published_relative_error = abs(computed - published) ./ published;
assert(max(published_relative_error) <= 0.005);

% Ritter solution: instantaneous release into a dry, horizontal rectangle.
g = 9.81;
width_m = 10;
reservoir_depth_m = 10;
dam_x_m = 1000;
final_time_s = 20;
cell_count = 400;
face_x_m = linspace(0, 2000, cell_count + 1).';
dx_m = face_x_m(2) - face_x_m(1);
cell_x_m = face_x_m(1:end-1) + dx_m / 2;
[cells, faces] = rectangular_grid(cell_x_m, face_x_m, width_m, 20);
area = zeros(cell_count, 1);
area(cell_x_m < dam_x_m) = width_m * reservoir_depth_m;
discharge = zeros(cell_count, 1);
initial_volume_m3 = sum(area) * dx_m;
dt_s = 0.30 * dx_m / sqrt(g * reservoir_depth_m);
time_s = 0;
total_adjustment_m3 = 0;
while time_s < final_time_s
    step_s = min(dt_s, final_time_s - time_s);
    [area, discharge, diagnostic] = hydrohp_varying_step( ...
        cells, faces, area, discharge, dx_m, step_s, ...
        include_friction=false, cfl_target=0.7, minimum_area_m2=1e-10);
    total_adjustment_m3 = total_adjustment_m3 + diagnostic.volume_adjustment_m3;
    time_s = time_s + step_s;
end

depth_m = area / width_m;
celerity_m_s = sqrt(g * reservoir_depth_m);
xi = (cell_x_m - dam_x_m) / final_time_s;
analytical_depth_m = zeros(cell_count, 1);
analytical_depth_m(xi <= -celerity_m_s) = reservoir_depth_m;
fan = xi > -celerity_m_s & xi < 2 * celerity_m_s;
analytical_depth_m(fan) = (2 * celerity_m_s - xi(fan)).^2 / (9 * g);
relative_l1_depth = sum(abs(depth_m - analytical_depth_m)) / ...
    sum(analytical_depth_m);
front_thresholds_m = [0.01, 0.05, 0.10];
front_error_cells_by_threshold = zeros(size(front_thresholds_m));
for threshold_index = 1:numel(front_thresholds_m)
    threshold_m = front_thresholds_m(threshold_index);
    numerical_wet = find(depth_m > threshold_m, 1, 'last');
    analytical_wet = find(analytical_depth_m > threshold_m, 1, 'last');
    front_error_cells_by_threshold(threshold_index) = ...
        abs(numerical_wet - analytical_wet);
end
front_error_cells = front_error_cells_by_threshold(2);
analytical_front_travel_m = 2 * celerity_m_s * final_time_s;
relative_front_error = front_error_cells * dx_m / analytical_front_travel_m;
final_volume_m3 = sum(area) * dx_m;
relative_mass_residual = (final_volume_m3 - initial_volume_m3 - ...
    total_adjustment_m3) / initial_volume_m3;

fprintf(['Ritter diagnostics: L1 %.12g, front errors [%.3g %.3g %.3g] ' ...
    'cells at [0.01 0.05 0.10] m, mass %.12g\n'], relative_l1_depth, ...
    front_error_cells_by_threshold, relative_mass_residual);
assert(relative_l1_depth <= 0.08);
assert(relative_front_error <= 0.15);
assert(abs(relative_mass_residual) <= 1e-10);
assert(all(isfinite([area; discharge])) && all(area >= 0));

metrics = struct( ...
    'published_triangle_discharge_m3s', discharge_m3s, ...
    'published_triangle_max_relative_error', max(published_relative_error), ...
    'ritter_relative_l1_depth', relative_l1_depth, ...
    'ritter_front_error_cells', front_error_cells, ...
    'ritter_relative_front_error', relative_front_error, ...
    'ritter_relative_mass_residual', relative_mass_residual);
fprintf(['HydroHP H7 checks passed: published max error %.3g%%; ' ...
    'Ritter L1 %.4g, front %.3g cells, mass %.3g.\n'], ...
    100 * max(published_relative_error), relative_l1_depth, ...
    front_error_cells, relative_mass_residual);
end

function [cells, faces] = rectangular_grid(cell_x_m, face_x_m, width_m, maximum_depth_m)
template = rectangular_section(width_m, maximum_depth_m);
cells = repmat(template, numel(cell_x_m), 1);
faces = repmat(template, numel(face_x_m), 1);
for index = 1:numel(cell_x_m)
    cells(index).x_m = cell_x_m(index);
end
for index = 1:numel(face_x_m)
    faces(index).x_m = face_x_m(index);
end
end

function section = rectangular_section(width_m, maximum_depth_m)
depth_m = linspace(0, maximum_depth_m, 401).';
area_m2 = width_m * depth_m;
wetted_perimeter_m = width_m + 2 * depth_m;
wetted_perimeter_m(1) = 0;
top_width_m = width_m * ones(size(depth_m));
top_width_m(1) = 0;
hydraulic_radius_m = zeros(size(depth_m));
wet = depth_m > 0;
hydraulic_radius_m(wet) = area_m2(wet) ./ wetted_perimeter_m(wet);
hydrostatic_integral_m3 = 0.5 * width_m * depth_m.^2;
centroid_depth_m = 0.5 * depth_m;
effective_manning_n = 0.03 * ones(size(depth_m));
conveyance_m3s = zeros(size(depth_m));
conveyance_m3s(wet) = area_m2(wet) .* hydraulic_radius_m(wet).^(2/3) / 0.03;
section = struct('x_m', 0, 'invert_elevation_m', 0, ...
    'maximum_water_surface_elevation_m', maximum_depth_m, ...
    'table', table(depth_m, depth_m, area_m2, wetted_perimeter_m, ...
    top_width_m, hydraulic_radius_m, hydrostatic_integral_m3, ...
    centroid_depth_m, effective_manning_n, conveyance_m3s, ...
    'VariableNames', {'depth_m', 'water_surface_elevation_m', 'area_m2', ...
    'wetted_perimeter_m', 'top_width_m', 'hydraulic_radius_m', ...
    'hydrostatic_integral_m3', 'centroid_depth_m', ...
    'effective_manning_n', 'conveyance_m3s'}));
end
