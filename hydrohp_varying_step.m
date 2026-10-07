function [area_new, discharge_new, diagnostics] = hydrohp_varying_step( ...
        cell_sections, face_sections, area_m2, discharge_m3s, dx_m, dt_s, options)
%HYDROHP_VARYING_STEP Advance the varying-section Saint-Venant state once.
%
% The finite-volume step reconstructs both neighbouring states on the same
% interface geometry. A matching hydrostatic-pressure source preserves a
% horizontal, motionless water surface over varying bed and section shape.

arguments
    cell_sections (:,1) struct
    face_sections (:,1) struct
    area_m2 (:,1) double
    discharge_m3s (:,1) double
    dx_m (1,1) double {mustBeFinite, mustBePositive}
    dt_s (1,1) double {mustBeFinite, mustBePositive}
    options.gravity_m_s2 (1,1) double {mustBeFinite, mustBePositive} = 9.81
    options.cfl_target (1,1) double {mustBeFinite, mustBePositive} = 0.8
    options.minimum_area_m2 (1,1) double {mustBeFinite, mustBeNonnegative} = 1e-10
    options.include_friction (1,1) logical = true
end

cell_count = numel(cell_sections);
if numel(face_sections) ~= cell_count + 1 || ...
        numel(area_m2) ~= cell_count || numel(discharge_m3s) ~= cell_count
    error('hydrohp:GridSizeMismatch', ...
        'The grid requires N cells, N states, and N+1 face sections.');
end
if any(~isfinite(area_m2)) || any(area_m2 < 0) || any(~isfinite(discharge_m3s))
    error('hydrohp:InvalidHydrodynamicState', ...
        'Area and discharge must be finite, and area must be nonnegative.');
end
if options.cfl_target > 1
    error('hydrohp:InvalidCFLTarget', 'cfl_target must not exceed 1.');
end

g = options.gravity_m_s2;
minimum_area = options.minimum_area_m2;
stage = zeros(cell_count, 1);
velocity = zeros(cell_count, 1);
cell_conveyance = zeros(cell_count, 1);
for cell_index = 1:cell_count
    state = state_from_area(cell_sections(cell_index), area_m2(cell_index));
    stage(cell_index) = state.stage_m;
    cell_conveyance(cell_index) = state.conveyance_m3s;
    if area_m2(cell_index) > minimum_area
        velocity(cell_index) = discharge_m3s(cell_index) / area_m2(cell_index);
    end
end

flux = zeros(cell_count + 1, 2);
pressure_at_right_face = zeros(cell_count, 1);
pressure_at_left_face = zeros(cell_count, 1);
maximum_face_speed = 0;

for face_index = 1:cell_count + 1
    if face_index == 1
        left_cell = 1;
        right_cell = 1;
        left_sign = -1;
        right_sign = 1;
    elseif face_index == cell_count + 1
        left_cell = cell_count;
        right_cell = cell_count;
        left_sign = 1;
        right_sign = -1;
    else
        left_cell = face_index - 1;
        right_cell = face_index;
        left_sign = 1;
        right_sign = 1;
    end

    left = state_from_stage(face_sections(face_index), stage(left_cell));
    right = state_from_stage(face_sections(face_index), stage(right_cell));
    q_left = left_sign * velocity(left_cell) * left.area_m2;
    q_right = right_sign * velocity(right_cell) * right.area_m2;
    [flux(face_index, :), face_speed] = rusanov_flux( ...
        left, right, q_left, q_right, g, minimum_area);
    maximum_face_speed = max(maximum_face_speed, face_speed);

    if face_index > 1
        pressure_at_right_face(face_index - 1) = left.hydrostatic_integral_m3;
    end
    if face_index <= cell_count
        pressure_at_left_face(face_index) = right.hydrostatic_integral_m3;
    end
end

if maximum_face_speed > 0
    cfl_number = dt_s * maximum_face_speed / dx_m;
    recommended_dt_s = options.cfl_target * dx_m / maximum_face_speed;
else
    cfl_number = 0;
    recommended_dt_s = inf;
end
if cfl_number > options.cfl_target * (1 + 10*eps)
    error('hydrohp:CFLViolation', ...
        'dt = %.6g s gives CFL %.6g; use dt <= %.6g s.', ...
        dt_s, cfl_number, recommended_dt_s);
end

area_raw = area_m2 - dt_s / dx_m * diff(flux(:, 1));
geometric_source = g / dx_m * ...
    (pressure_at_right_face - pressure_at_left_face);
discharge_raw = discharge_m3s - dt_s / dx_m * diff(flux(:, 2)) + ...
    dt_s * geometric_source;

if options.include_friction
    wet = area_raw > minimum_area & cell_conveyance > 0;
    friction_coefficient = zeros(cell_count, 1);
    friction_coefficient(wet) = g * area_raw(wet) ./ cell_conveyance(wet).^2;
    discharge_raw(wet) = discharge_raw(wet) ./ ...
        (1 + dt_s * friction_coefficient(wet) .* abs(discharge_raw(wet)));
end

area_new = area_raw;
discharge_new = discharge_raw;
dry = area_new <= minimum_area;
area_new(dry) = 0;
discharge_new(dry) = 0;
volume_adjustment_m3 = sum(area_new - area_raw) * dx_m;

if any(~isfinite(area_new)) || any(~isfinite(discharge_new)) || any(area_new < 0)
    error('hydrohp:InvalidUpdatedState', ...
        'The varying-section update produced a nonfinite or negative state.');
end

diagnostics = struct( ...
    'cfl_number', cfl_number, ...
    'maximum_wave_speed_m_s', maximum_face_speed, ...
    'recommended_dt_s', recommended_dt_s, ...
    'boundary_inflow_volume_m3', dt_s * flux(1, 1), ...
    'boundary_outflow_volume_m3', dt_s * flux(end, 1), ...
    'volume_adjustment_m3', volume_adjustment_m3, ...
    'storage_change_m3', sum(area_new - area_m2) * dx_m, ...
    'geometric_source_m4_s2', geometric_source);
end

function state = state_from_area(section, area_m2)
table_data = validate_section(section);
maximum_area = table_data.area_m2(end);
area_tolerance = 1e-10 * max(1, maximum_area);
if area_m2 > maximum_area + area_tolerance
    error('hydrohp:AreaOutsideSectionRange', ...
        'Area %.6g m2 exceeds the local section range.', area_m2);
end
area_m2 = min(area_m2, maximum_area);
depth_m = interp1(table_data.area_m2, table_data.depth_m, area_m2, 'pchip');
state = state_from_stage(section, section.invert_elevation_m + depth_m);
end

function state = state_from_stage(section, stage_m)
table_data = validate_section(section);
if stage_m <= section.invert_elevation_m
    state = struct('stage_m', section.invert_elevation_m, 'area_m2', 0, ...
        'top_width_m', 0, 'hydrostatic_integral_m3', 0, 'conveyance_m3s', 0);
    return
end
stage_tolerance = 1e-10 * max(1, ...
    abs(section.maximum_water_surface_elevation_m));
if stage_m > section.maximum_water_surface_elevation_m + stage_tolerance
    error('hydrohp:StageOutsideSectionRange', ...
        'Stage %.6g m exceeds the local section range.', stage_m);
end
stage_m = min(stage_m, section.maximum_water_surface_elevation_m);
state = struct( ...
    'stage_m', stage_m, ...
    'area_m2', interp1(table_data.water_surface_elevation_m, ...
        table_data.area_m2, stage_m, 'pchip'), ...
    'top_width_m', interp1(table_data.water_surface_elevation_m, ...
        table_data.top_width_m, stage_m, 'linear'), ...
    'hydrostatic_integral_m3', interp1(table_data.water_surface_elevation_m, ...
        table_data.hydrostatic_integral_m3, stage_m, 'pchip'), ...
    'conveyance_m3s', interp1(table_data.water_surface_elevation_m, ...
        table_data.conveyance_m3s, stage_m, 'pchip'));
end

function table_data = validate_section(section)
required = {'invert_elevation_m', 'maximum_water_surface_elevation_m', 'table'};
table_fields = {'depth_m', 'water_surface_elevation_m', 'area_m2', ...
    'top_width_m', 'hydrostatic_integral_m3', 'conveyance_m3s'};
if ~isscalar(section) || ~all(isfield(section, required)) || ...
        ~istable(section.table) || ...
        ~all(ismember(table_fields, section.table.Properties.VariableNames))
    error('hydrohp:InvalidSectionTable', ...
        'Each local section must contain an H4/H5a hydraulic table.');
end
table_data = section.table;
end

function [flux, maximum_speed] = rusanov_flux( ...
        left, right, q_left, q_right, g, minimum_area)
u_left = 0;
u_right = 0;
c_left = 0;
c_right = 0;
if left.area_m2 > minimum_area
    u_left = q_left / left.area_m2;
    c_left = sqrt(g * left.area_m2 / max(left.top_width_m, eps));
end
if right.area_m2 > minimum_area
    u_right = q_right / right.area_m2;
    c_right = sqrt(g * right.area_m2 / max(right.top_width_m, eps));
end
maximum_speed = max(abs(u_left) + c_left, abs(u_right) + c_right);
physical_left = [q_left, q_left * u_left + g * left.hydrostatic_integral_m3];
physical_right = [q_right, q_right * u_right + g * right.hydrostatic_integral_m3];
flux = 0.5 * (physical_left + physical_right) - 0.5 * maximum_speed * ...
    ([right.area_m2, q_right] - [left.area_m2, q_left]);
end
