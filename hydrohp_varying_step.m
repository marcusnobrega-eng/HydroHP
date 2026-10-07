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
    options.boundary (1,1) struct = struct()
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
        inner = state_from_stage(face_sections(face_index), stage(1));
        inner_discharge = velocity(1) * inner.area_m2;
        pressure_at_left_face(1) = inner.hydrostatic_integral_m3;
        if isfield(options.boundary, 'upstream_discharge_m3s')
            [flux(face_index, :), face_speed] = upstream_flux( ...
                face_sections(face_index), inner, inner_discharge, ...
                options.boundary.upstream_discharge_m3s, g, minimum_area);
        else
            [flux(face_index, :), face_speed] = rusanov_flux( ...
                inner, inner, -inner_discharge, inner_discharge, g, minimum_area);
        end
    elseif face_index == cell_count + 1
        inner = state_from_stage(face_sections(face_index), stage(cell_count));
        inner_discharge = velocity(cell_count) * inner.area_m2;
        pressure_at_right_face(cell_count) = inner.hydrostatic_integral_m3;
        if isfield(options.boundary, 'downstream_type')
            [flux(face_index, :), face_speed] = downstream_flux( ...
                face_sections(face_index), inner, inner_discharge, ...
                options.boundary, g, minimum_area);
        else
            [flux(face_index, :), face_speed] = rusanov_flux( ...
                inner, inner, inner_discharge, -inner_discharge, g, minimum_area);
        end
    else
        left_cell = face_index - 1;
        right_cell = face_index;
        left = state_from_stage(face_sections(face_index), stage(left_cell));
        right = state_from_stage(face_sections(face_index), stage(right_cell));
        q_left = velocity(left_cell) * left.area_m2;
        q_right = velocity(right_cell) * right.area_m2;
        [flux(face_index, :), face_speed] = rusanov_flux( ...
            left, right, q_left, q_right, g, minimum_area);
        pressure_at_right_face(face_index - 1) = left.hydrostatic_integral_m3;
        pressure_at_left_face(face_index) = right.hydrostatic_integral_m3;
    end
    maximum_face_speed = max(maximum_face_speed, face_speed);
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
        (1 + dt_s * friction_coefficient(wet) .* abs(discharge_m3s(wet)));
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

function [flux, speed] = upstream_flux( ...
        section, inner, inner_discharge, prescribed_discharge, g, minimum_area)
if ~isscalar(prescribed_discharge) || ~isfinite(prescribed_discharge) || ...
        prescribed_discharge < 0
    error('hydrohp:InvalidUpstreamDischarge', ...
        'The prescribed upstream discharge must be a finite nonnegative scalar.');
end
[inner_phi, inner_celerity] = characteristic_values(section, inner.area_m2, g);
inner_velocity = safe_velocity(inner_discharge, inner.area_m2, minimum_area);
if inner_velocity - inner_celerity >= 0 && inner.area_m2 > minimum_area
    error('hydrohp:UnsupportedBoundaryRegime', ...
        'Discharge alone is insufficient at a supercritical upstream inflow.');
end
invariant_minus = inner_velocity - inner_phi;
residual = @(area) prescribed_discharge ./ area - ...
    characteristic_values(section, area, g) - invariant_minus;
state = solve_subcritical_state( ...
    section, residual, @(area) prescribed_discharge + zeros(size(area)), ...
    inner.area_m2, g, minimum_area);
flux = physical_flux(state, prescribed_discharge, g, minimum_area);
speed = abs(safe_velocity(prescribed_discharge, state.area_m2, minimum_area)) + ...
    wave_celerity(state, g, minimum_area);
end

function [flux, speed] = downstream_flux( ...
        section, inner, inner_discharge, boundary, g, minimum_area)
[inner_phi, inner_celerity] = characteristic_values(section, inner.area_m2, g);
inner_velocity = safe_velocity(inner_discharge, inner.area_m2, minimum_area);
if inner_velocity - inner_celerity >= 0
    flux = physical_flux(inner, inner_discharge, g, minimum_area);
    speed = abs(inner_velocity) + inner_celerity;
    return
end
if inner_velocity + inner_celerity <= 0 && inner.area_m2 > minimum_area
    error('hydrohp:UnsupportedBoundaryRegime', ...
        'A downstream inflow requires more than one boundary condition.');
end
invariant_plus = inner_velocity + inner_phi;
boundary_type = lower(string(boundary.downstream_type));
if boundary_type == "normal_depth"
    if ~isfield(boundary, 'downstream_slope') || ...
            ~isscalar(boundary.downstream_slope) || ...
            ~isfinite(boundary.downstream_slope) || boundary.downstream_slope <= 0
        error('hydrohp:InvalidDownstreamSlope', ...
            'A normal-depth boundary requires a positive downstream_slope.');
    end
    discharge = @(area) conveyance_at_area(section, area) * ...
        sqrt(boundary.downstream_slope);
    residual = @(area) discharge(area) ./ area + ...
        characteristic_values(section, area, g) - invariant_plus;
    state = solve_subcritical_state( ...
        section, residual, discharge, inner.area_m2, g, minimum_area);
    boundary_discharge = discharge(state.area_m2);
elseif boundary_type == "stage"
    if ~isfield(boundary, 'downstream_stage_m') || ...
            ~isscalar(boundary.downstream_stage_m) || ...
            ~isfinite(boundary.downstream_stage_m)
        error('hydrohp:InvalidDownstreamStage', ...
            'A stage boundary requires downstream_stage_m.');
    end
    state = state_from_stage(section, boundary.downstream_stage_m);
    phi = characteristic_values(section, state.area_m2, g);
    boundary_discharge = state.area_m2 * (invariant_plus - phi);
    if abs(safe_velocity(boundary_discharge, state.area_m2, minimum_area)) >= ...
            wave_celerity(state, g, minimum_area)
        error('hydrohp:UnsupportedBoundaryRegime', ...
            'The prescribed downstream stage creates a non-subcritical state.');
    end
else
    error('hydrohp:InvalidDownstreamType', ...
        'downstream_type must be "normal_depth" or "stage".');
end
flux = physical_flux(state, boundary_discharge, g, minimum_area);
speed = abs(safe_velocity(boundary_discharge, state.area_m2, minimum_area)) + ...
    wave_celerity(state, g, minimum_area);
end

function state = solve_subcritical_state( ...
        section, residual, discharge, preferred_area, g, minimum_area)
areas = section.table.area_m2;
wet = areas > minimum_area;
areas = areas(wet);
values = residual(areas);
candidates = zeros(0, 1);
for index = 1:numel(areas) - 1
    if values(index) == 0
        candidates(end + 1, 1) = areas(index); %#ok<AGROW>
    elseif isfinite(values(index)) && isfinite(values(index + 1)) && ...
            sign(values(index)) ~= sign(values(index + 1))
        candidates(end + 1, 1) = fzero(residual, ...
            [areas(index), areas(index + 1)]); %#ok<AGROW>
    end
end
if values(end) == 0
    candidates(end + 1, 1) = areas(end);
end
if ~isempty(candidates)
    q = discharge(candidates);
    states = arrayfun(@(area) state_from_area(section, area), candidates);
    celerity = arrayfun(@(item) wave_celerity(item, g, minimum_area), states);
    froude = abs(q ./ candidates) ./ celerity;
    candidates = candidates(froude < 1);
end
if isempty(candidates)
    error('hydrohp:BoundaryStateUnavailable', ...
        'No subcritical boundary state exists within the section table.');
end
[~, selected] = min(abs(candidates - preferred_area));
state = state_from_area(section, candidates(selected));
end

function [phi, celerity] = characteristic_values(section, area, g)
table_data = section.table;
wet = table_data.area_m2 > 0 & table_data.top_width_m > 0;
depth = table_data.depth_m(wet);
table_area = table_data.area_m2(wet);
integrand = sqrt(g * table_data.top_width_m(wet) ./ table_area);
table_phi = zeros(size(table_area));
table_phi(1) = 2 * depth(1) * integrand(1);
if numel(depth) > 1
    integrated = cumtrapz(depth, integrand);
    table_phi(2:end) = table_phi(1) + integrated(2:end);
end
phi = zeros(size(area));
celerity = zeros(size(area));
positive = area > 0;
if any(positive)
    phi(positive) = interp1(table_area, table_phi, area(positive), 'pchip');
    width = interp1(table_data.area_m2, table_data.top_width_m, ...
        area(positive), 'linear');
    celerity(positive) = sqrt(g * area(positive) ./ max(width, eps));
end
end

function conveyance = conveyance_at_area(section, area)
conveyance = interp1(section.table.area_m2, ...
    section.table.conveyance_m3s, area, 'pchip');
end

function flux = physical_flux(state, discharge, g, minimum_area)
velocity = safe_velocity(discharge, state.area_m2, minimum_area);
flux = [discharge, discharge * velocity + ...
    g * state.hydrostatic_integral_m3];
end

function velocity = safe_velocity(discharge, area, minimum_area)
velocity = 0;
if area > minimum_area
    velocity = discharge / area;
end
end

function celerity = wave_celerity(state, g, minimum_area)
celerity = 0;
if state.area_m2 > minimum_area
    celerity = sqrt(g * state.area_m2 / max(state.top_width_m, eps));
end
end
