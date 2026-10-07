function test_varying_geometry_kernel
%TEST_VARYING_GEOMETRY_KERNEL Check H5b balance, CFL, wet/dry, and volume.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);

% Lake at rest over a raised bed and smooth width contraction.
[cells, faces, dx_m] = make_grid([10; 6; 10], [0; 0.4; 0], 40);
stage_m = 1.5 * ones(numel(cells), 1);
area = areas_at_stage(cells, stage_m);
discharge = zeros(size(area));
[area_next, discharge_next, diagnostic] = hydrohp_varying_step( ...
    cells, faces, area, discharge, dx_m, 0.05, ...
    include_friction=false, cfl_target=0.8);
assert(max(abs(area_next - area)) < 1e-13);
assert(max(abs(discharge_next)) < 1e-12);
assert(abs(diagnostic.storage_change_m3) < 1e-12);

% A small wave must cross the contraction without losing closed-domain volume.
x = [cells.x_m].';
stage_m = 1.0 + 0.08 * exp(-((x - 25) / 8).^2);
area = areas_at_stage(cells, stage_m);
discharge = zeros(size(area));
initial_volume = sum(area) * dx_m;
total_adjustment = 0;
for step = 1:120
    dt_s = stable_dt(cells, faces, area, discharge, dx_m);
    [area, discharge, diagnostic] = hydrohp_varying_step( ...
        cells, faces, area, discharge, dx_m, dt_s, ...
        include_friction=false, cfl_target=0.7);
    total_adjustment = total_adjustment + diagnostic.volume_adjustment_m3;
end
final_volume = sum(area) * dx_m;
assert(all(isfinite([area; discharge])) && all(area >= 0));
assert(abs((final_volume - initial_volume) - total_adjustment) < 1e-10);

% A dam-break front into a dry bed must remain nonnegative and conserve volume.
[cells, faces, dx_m] = make_grid([8; 8; 8], [0; 0; 0], 60);
stage_m = zeros(numel(cells), 1);
stage_m([cells.x_m].' < 50) = 1;
area = areas_at_stage(cells, stage_m);
discharge = zeros(size(area));
initial_volume = sum(area) * dx_m;
total_adjustment = 0;
for step = 1:80
    dt_s = stable_dt(cells, faces, area, discharge, dx_m);
    [area, discharge, diagnostic] = hydrohp_varying_step( ...
        cells, faces, area, discharge, dx_m, dt_s, ...
        include_friction=false, cfl_target=0.7, minimum_area_m2=1e-9);
    total_adjustment = total_adjustment + diagnostic.volume_adjustment_m3;
end
final_volume = sum(area) * dx_m;
assert(all(isfinite([area; discharge])) && all(area >= 0));
assert(any(area([cells.x_m].' >= 50) > 0));
assert(abs((final_volume - initial_volume) - total_adjustment) < 1e-9);

% The local conveyance table must damp a moving wet state when friction is on.
uniform_stage = ones(numel(cells), 1);
uniform_area = areas_at_stage(cells, uniform_stage);
uniform_discharge = 2 * ones(size(uniform_area));
[~, no_friction] = hydrohp_varying_step( ...
    cells, faces, uniform_area, uniform_discharge, dx_m, 0.01, ...
    include_friction=false, cfl_target=0.7);
[~, with_friction] = hydrohp_varying_step( ...
    cells, faces, uniform_area, uniform_discharge, dx_m, 0.01, ...
    include_friction=true, cfl_target=0.7);
assert(abs(with_friction(10)) < abs(no_friction(10)));

% Unsafe time steps must fail explicitly.
try
    hydrohp_varying_step(cells, faces, area, discharge, dx_m, 100, ...
        include_friction=false, cfl_target=0.7);
    error('hydrohp:TestFailed', 'An unsafe CFL step was accepted.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:CFLViolation'));
end

fprintf('HydroHP varying-geometry kernel checks passed.\n');
end

function [cells, faces, dx_m] = make_grid(width_m, invert_m, cell_count)
chainage_m = [0; 50; 100];
rows = table();
for section_index = 1:numel(chainage_m)
    section_rows = table( ...
        repmat("S" + section_index, 4, 1), ...
        repmat(chainage_m(section_index), 4, 1), ...
        width_m(section_index) * [0; 0.1; 0.9; 1], ...
        invert_m(section_index) + [3; 0; 0; 3], ...
        0.03 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section_rows]; %#ok<AGROW>
end
surveyed = hydrohp_build_section_tables(rows, 0.02, "scm");
face_x = linspace(0, 100, cell_count + 1).';
dx_m = face_x(2) - face_x(1);
cell_x = face_x(1:end-1) + dx_m / 2;
cells = hydrohp_map_sections_to_grid(surveyed, cell_x);
faces = hydrohp_map_sections_to_grid(surveyed, face_x);
end

function area = areas_at_stage(sections, stage_m)
area = zeros(numel(sections), 1);
for index = 1:numel(sections)
    if stage_m(index) > sections(index).invert_elevation_m
        area(index) = interp1( ...
            sections(index).table.water_surface_elevation_m, ...
            sections(index).table.area_m2, stage_m(index), 'pchip');
    end
end
end

function dt_s = stable_dt(cells, faces, area, discharge, dx_m)
% Probe with a very small step; the kernel returns its current CFL limit.
[~, ~, diagnostic] = hydrohp_varying_step( ...
    cells, faces, area, discharge, dx_m, 1e-6, ...
    include_friction=false, cfl_target=0.7);
dt_s = min(0.05, diagnostic.recommended_dt_s);
end
