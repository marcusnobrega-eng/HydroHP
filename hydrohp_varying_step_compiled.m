function [area_new, discharge_new, diagnostics] = hydrohp_varying_step_compiled( ...
        geometry, area_m2, discharge_m3s, dx_m, dt_s, options)
%HYDROHP_VARYING_STEP_COMPILED Vectorized form of HYDROHP_VARYING_STEP.
%
% The update is the same finite-volume step as HYDROHP_VARYING_STEP, applied
% to geometry from HYDROHP_COMPILE_GEOMETRY. Cell and interior-face states
% are evaluated for all locations at once; the two boundary faces use the
% scalar boundary solvers of the reference kernel. The reference kernel is
% kept unchanged so tests can require identical results.
%
% reconstruction = "stage" (default) reconstructs both neighbours of an
% interior face from their water-surface elevation on the face section,
% exactly as the reference kernel. reconstruction = "hydrostatic" applies
% the first-order hydrostatic reconstruction of Audusse et al. (2004): the
% face bed is max(z_L, z_R) of the neighbouring cell inverts and each
% reconstructed depth is max(0, eta - max(z_L, z_R)), evaluated on the face
% section. The reconstructed depth never exceeds the cell depth, which
% restores the positivity argument of the HLL update, and equal surfaces
% still give equal face states, so lake at rest is preserved. The scheme is
% known to overestimate depth where the bed drop per cell exceeds the water
% depth (Delestre et al., 2012); those faces are counted in the diagnostics.

arguments
    geometry (1,1) struct
    area_m2 (:,1) double
    discharge_m3s (:,1) double
    dx_m (1,1) double {mustBeFinite, mustBePositive}
    dt_s (1,1) double {mustBeFinite, mustBePositive}
    options.gravity_m_s2 (1,1) double {mustBeFinite, mustBePositive} = 9.81
    options.cfl_target (1,1) double {mustBeFinite, mustBePositive} = 0.8
    options.minimum_area_m2 (1,1) double {mustBeFinite, mustBeNonnegative} = 1e-10
    options.include_friction (1,1) logical = true
    options.boundary (1,1) struct = struct()
    options.report_negative_area (1,1) logical = false
    options.negative_area_tolerance_m2 (1,1) double {mustBeFinite, mustBeNonnegative} = 1e-10
    options.reconstruction (1,1) string {mustBeMember(options.reconstruction, ...
        ["stage", "hydrostatic"])} = "stage"
end

cell_count = geometry.cell_count;
if numel(area_m2) ~= cell_count || numel(discharge_m3s) ~= cell_count
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
if g ~= geometry.gravity_m_s2
    error('hydrohp:GravityMismatch', ...
        'The geometry was compiled for g = %.6g m/s2, not %.6g m/s2.', ...
        geometry.gravity_m_s2, g);
end
minimum_area = options.minimum_area_m2;
cells = geometry.cells;
faces = geometry.faces;

[stage, cell_conveyance] = cell_states(cells, area_m2);
velocity = zeros(cell_count, 1);
wet_cell = area_m2 > minimum_area;
velocity(wet_cell) = discharge_m3s(wet_cell) ./ area_m2(wet_cell);

% Interior faces 2..N: both neighbours reconstructed on the face geometry.
interior = (2:cell_count).';
if options.reconstruction == "stage"
    left_query = stage(interior - 1);
    right_query = stage(interior);
else
    face_bed = max(cells.invert_elevation_m(interior - 1), ...
        cells.invert_elevation_m(interior));
    left_depth = max(stage(interior - 1) - face_bed, 0);
    right_depth = max(stage(interior) - face_bed, 0);
    face_invert = faces.invert_elevation_m(interior);
    left_query = face_invert + left_depth;
    right_query = face_invert + right_depth;
    cell_depth = stage - cells.invert_elevation_m;
    bed_drop = cells.invert_elevation_m(interior) - ...
        cells.invert_elevation_m(interior - 1);
    steep_shallow_face_count = sum( ...
        (bed_drop >= cell_depth(interior - 1) & cell_depth(interior - 1) > 0) | ...
        (-bed_drop >= cell_depth(interior) & cell_depth(interior) > 0));
end
[left, left_excess_m] = reconstruct(faces, interior, left_query);
[right, right_excess_m] = reconstruct(faces, interior, right_query);
q_left = velocity(interior - 1) .* left.area_m2;
q_right = velocity(interior) .* right.area_m2;
[interior_flux, interior_speed] = hll_flux_vector( ...
    left, right, q_left, q_right, g, minimum_area);

% Boundary faces use the scalar boundary solvers.
upstream_section = geometry.upstream_face;
[upstream_inner, upstream_excess_m] = boundary_inner_state( ...
    upstream_section, cells.invert_elevation_m(1), stage(1), ...
    options.boundary, "upstream");
upstream_inner_discharge = velocity(1) * upstream_inner.area_m2;
if isfield(options.boundary, 'upstream_discharge_m3s')
    [upstream_flux_value, upstream_speed] = upstream_flux( ...
        upstream_section, upstream_inner, upstream_inner_discharge, ...
        options.boundary.upstream_discharge_m3s, ...
        options.boundary, g, minimum_area);
else
    [upstream_flux_value, upstream_speed] = hll_flux_vector( ...
        upstream_inner, upstream_inner, -upstream_inner_discharge, ...
        upstream_inner_discharge, g, minimum_area);
end

downstream_section = geometry.downstream_face;
[downstream_inner, downstream_excess_m] = boundary_inner_state( ...
    downstream_section, cells.invert_elevation_m(cell_count), ...
    stage(cell_count), options.boundary, "downstream");
downstream_inner_discharge = velocity(cell_count) * downstream_inner.area_m2;
if isfield(options.boundary, 'downstream_type')
    [downstream_flux_value, downstream_speed] = downstream_flux( ...
        downstream_section, downstream_inner, downstream_inner_discharge, ...
        options.boundary, g, minimum_area);
else
    [downstream_flux_value, downstream_speed] = hll_flux_vector( ...
        downstream_inner, downstream_inner, downstream_inner_discharge, ...
        -downstream_inner_discharge, g, minimum_area);
end

flux = [upstream_flux_value; interior_flux; downstream_flux_value];
face_speed = [upstream_speed; interior_speed; downstream_speed];
pressure_at_right_face = [left.hydrostatic_integral_m3; ...
    downstream_inner.hydrostatic_integral_m3];
pressure_at_left_face = [upstream_inner.hydrostatic_integral_m3; ...
    right.hydrostatic_integral_m3];
all_excess_m = [upstream_excess_m; left_excess_m; right_excess_m; ...
    downstream_excess_m];
stage_endpoint_clamp_count = sum(all_excess_m > 0);
maximum_stage_endpoint_excess_m = max([0; all_excess_m]);

[maximum_face_speed, maximum_face_index] = max(face_speed);
if ~(maximum_face_speed > 0)
    maximum_face_speed = 0;
    maximum_face_index = 1;
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
        ['dt = %.6g s gives CFL %.6g at face %d (x = %.6g m); ' ...
        'use dt <= %.6g s.'], dt_s, cfl_number, maximum_face_index, ...
        faces.x_m(maximum_face_index), recommended_dt_s);
end

area_raw = area_m2 - dt_s / dx_m * diff(flux(:, 1));
geometric_source = g / dx_m * ...
    (pressure_at_right_face - pressure_at_left_face);
discharge_raw = discharge_m3s - dt_s / dx_m * diff(flux(:, 2)) + ...
    dt_s * geometric_source;
discharge_before_friction = discharge_raw;
friction_coefficient = zeros(cell_count, 1);

if options.include_friction
    wet = area_raw > minimum_area & cell_conveyance > 0;
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
    'maximum_wave_speed_face_index', maximum_face_index, ...
    'recommended_dt_s', recommended_dt_s, ...
    'boundary_inflow_volume_m3', dt_s * flux(1, 1), ...
    'boundary_outflow_volume_m3', dt_s * flux(end, 1), ...
    'volume_adjustment_m3', volume_adjustment_m3, ...
    'stage_endpoint_clamp_count', stage_endpoint_clamp_count, ...
    'maximum_stage_endpoint_excess_m', maximum_stage_endpoint_excess_m, ...
    'storage_change_m3', sum(area_new - area_m2) * dx_m, ...
    'geometric_source_m4_s2', geometric_source);
if options.reconstruction == "hydrostatic"
    diagnostics.steep_shallow_face_count = steep_shallow_face_count;
end

% Opt-in P0.1 diagnostic: describe raw areas below -tolerance before they are
% clipped. It reads the step's intermediate values and never changes them.
if options.report_negative_area
    event = struct([]);
    negative = find(area_raw < -options.negative_area_tolerance_m2);
    if ~isempty(negative)
        face_left = join_states(upstream_inner, left, downstream_inner);
        face_right = join_states(upstream_inner, right, downstream_inner);
        event = negative_area_event(negative, cells, faces, area_m2, ...
            discharge_m3s, stage, velocity, cell_conveyance, area_raw, ...
            discharge_before_friction, discharge_raw, friction_coefficient, ...
            flux, face_speed, face_left, face_right, pressure_at_left_face, ...
            pressure_at_right_face, geometric_source, dx_m, dt_s, g, ...
            diagnostics);
    end
    diagnostics.negative_area_event = event;
end
end

% ---------------------------------------------------------------------------
% Negative raw-area diagnostic (P0.1)

function states = join_states(upstream, interior, downstream)
names = {'stage_m', 'area_m2', 'top_width_m', 'hydrostatic_integral_m3'};
states = struct();
for name = names
    states.(name{1}) = [upstream.(name{1}); interior.(name{1}); ...
        downstream.(name{1})];
end
end

function event = negative_area_event(negative, cells, faces, area, ...
        discharge, stage, velocity, conveyance, area_raw, ...
        discharge_before_friction, discharge_raw, friction_coefficient, ...
        flux, face_speed, face_left, face_right, pressure_left, ...
        pressure_right, geometric_source, dx_m, dt_s, g, diagnostics)
[minimum_raw_area, most_negative] = min(area_raw(negative));
selected = unique([negative(1); negative(most_negative)], 'stable');
records = arrayfun(@(index) cell_record(index, cells, faces, area, ...
    discharge, stage, velocity, conveyance, area_raw, ...
    discharge_before_friction, discharge_raw, friction_coefficient, flux, ...
    face_speed, face_left, face_right, pressure_left, pressure_right, ...
    geometric_source, dx_m, dt_s, g), selected);
event = struct( ...
    'negative_cell_index', negative, ...
    'negative_cell_x_m', cells.x_m(negative), ...
    'negative_cell_count', numel(negative), ...
    'minimum_raw_area_m2', minimum_raw_area, ...
    'most_negative_cell_index', negative(most_negative), ...
    'negative_raw_volume_m3', -sum(area_raw(negative)) * dx_m, ...
    'dt_s', dt_s, ...
    'dx_m', dx_m, ...
    'cfl_number', diagnostics.cfl_number, ...
    'maximum_wave_speed_face_index', diagnostics.maximum_wave_speed_face_index, ...
    'boundary_inflow_volume_m3', diagnostics.boundary_inflow_volume_m3, ...
    'boundary_outflow_volume_m3', diagnostics.boundary_outflow_volume_m3, ...
    'volume_adjustment_m3', diagnostics.volume_adjustment_m3, ...
    'storage_before_m3', sum(area) * dx_m, ...
    'storage_after_raw_m3', sum(area_raw) * dx_m, ...
    'records', records);
end

function record = cell_record(index, cells, faces, area, discharge, ...
        stage, velocity, conveyance, area_raw, discharge_before_friction, ...
        discharge_raw, friction_coefficient, flux, face_speed, face_left, ...
        face_right, pressure_left, pressure_right, geometric_source, dx_m, ...
        dt_s, g)
cell_state = stage_states(cells, index, stage(index), true);
celerity = 0;
if cell_state.area_m2 > 0
    celerity = sqrt(g * cell_state.area_m2 / max(cell_state.top_width_m, eps));
end
neighbours = max(index - 1, 1):min(index + 1, cells.count);
left_face = index;
right_face = index + 1;
outgoing_volume = dt_s * (max(flux(right_face, 1), 0) + ...
    max(-flux(left_face, 1), 0));
record = struct( ...
    'cell_index', index, ...
    'x_m', cells.x_m(index), ...
    'invert_elevation_m', cells.invert_elevation_m(index), ...
    'previous_area_m2', area(index), ...
    'previous_discharge_m3s', discharge(index), ...
    'previous_stage_m', stage(index), ...
    'previous_depth_m', stage(index) - cells.invert_elevation_m(index), ...
    'previous_velocity_m_s', velocity(index), ...
    'previous_froude', abs(velocity(index)) / max(celerity, eps), ...
    'cell_top_width_m', cell_state.top_width_m, ...
    'cell_conveyance_m3s', conveyance(index), ...
    'area_raw_m2', area_raw(index), ...
    'discharge_before_friction_m3s', discharge_before_friction(index), ...
    'discharge_raw_m3s', discharge_raw(index), ...
    'friction_coefficient', friction_coefficient(index), ...
    'mass_flux_left_face_m3s', flux(left_face, 1), ...
    'mass_flux_right_face_m3s', flux(right_face, 1), ...
    'momentum_flux_left_face', flux(left_face, 2), ...
    'momentum_flux_right_face', flux(right_face, 2), ...
    'pressure_left_face_m3', pressure_left(index), ...
    'pressure_right_face_m3', pressure_right(index), ...
    'geometric_source_m4_s2', geometric_source(index), ...
    'face_speed_left_m_s', face_speed(left_face), ...
    'face_speed_right_m_s', face_speed(right_face), ...
    'local_cfl_left', dt_s * face_speed(left_face) / dx_m, ...
    'local_cfl_right', dt_s * face_speed(right_face) / dx_m, ...
    'outgoing_volume_m3', outgoing_volume, ...
    'available_volume_m3', area(index) * dx_m, ...
    'outgoing_to_available_ratio', outgoing_volume / max(area(index) * dx_m, eps), ...
    'left_face', face_record(faces, left_face, face_left, face_right), ...
    'right_face', face_record(faces, right_face, face_left, face_right), ...
    'neighbour_index', neighbours(:), ...
    'neighbour_area_m2', area(neighbours), ...
    'neighbour_discharge_m3s', discharge(neighbours), ...
    'neighbour_stage_m', stage(neighbours), ...
    'neighbour_invert_elevation_m', cells.invert_elevation_m(neighbours));
end

function record = face_record(faces, index, face_left, face_right)
record = struct( ...
    'face_index', index, ...
    'x_m', faces.x_m(index), ...
    'invert_elevation_m', faces.invert_elevation_m(index), ...
    'maximum_water_surface_elevation_m', ...
        faces.maximum_water_surface_elevation_m(index), ...
    'left_stage_m', face_left.stage_m(index), ...
    'left_area_m2', face_left.area_m2(index), ...
    'left_top_width_m', face_left.top_width_m(index), ...
    'left_hydrostatic_integral_m3', face_left.hydrostatic_integral_m3(index), ...
    'right_stage_m', face_right.stage_m(index), ...
    'right_area_m2', face_right.area_m2(index), ...
    'right_top_width_m', face_right.top_width_m(index), ...
    'right_hydrostatic_integral_m3', face_right.hydrostatic_integral_m3(index));
end

% ---------------------------------------------------------------------------
% Vectorized table evaluation

function [stage, conveyance] = cell_states(cells, area_m2)
% Cell stage and conveyance from area (state_from_area, then state_from_stage).
index = (1:cells.count).';
maximum_area = cells.maximum_area_m2;
area_tolerance = 1e-10 * max(1, maximum_area);
outside = area_m2 > maximum_area + area_tolerance;
if any(outside)
    first = find(outside, 1);
    error('hydrohp:AreaOutsideSectionRange', ...
        'Area %.6g m2 exceeds the local section range.', area_m2(first));
end
area_m2 = min(area_m2, maximum_area);
interval = find_interval(cells.area_m2, cells.rows, index, area_m2);
depth_m = horner(cells.depth_of_area, cells, index, interval, ...
    area_m2 - cells.area_m2(sub2ind(size(cells.area_m2), index, interval)));
state = stage_states(cells, index, cells.invert_elevation_m + depth_m, false);
stage = state.stage_m;
conveyance = state.conveyance_m3s;
end

function [state, stage_excess_m] = reconstruct(set, index, stage_m)
% reconstructed_state_from_stage for many face/stage pairs.
maximum_wse = set.maximum_water_surface_elevation_m(index);
stage_excess_m = max(stage_m - maximum_wse, 0);
beyond = stage_excess_m > set.table_resolution_m(index);
query = stage_m;
query(~beyond) = min(stage_m(~beyond), maximum_wse(~beyond));
state = stage_states(set, index, query, true);
end

function state = stage_states(set, index, stage_m, face_properties)
% state_from_stage for many section/stage pairs.
invert = set.invert_elevation_m(index);
maximum_wse = set.maximum_water_surface_elevation_m(index);
dry = stage_m <= invert;
stage_tolerance = 1e-10 * max(1, abs(maximum_wse));
outside = ~dry & stage_m > maximum_wse + stage_tolerance;
if any(outside)
    first = find(outside, 1);
    error('hydrohp:StageOutsideSectionRange', ...
        'Stage %.6g m exceeds the local section range [%.6g, %.6g] m at x = %.6g m.', ...
        stage_m(first), invert(first), maximum_wse(first), ...
        set.x_m(index(first)));
end
stage_m = min(stage_m, maximum_wse);
stage_m(dry) = invert(dry);
interval = find_interval(set.stage_m, set.rows, index, stage_m);
breakpoint = sub2ind(size(set.stage_m), index, interval);
local = stage_m - set.stage_m(breakpoint);
conveyance = horner(set.conveyance_of_stage, set, index, interval, local);
conveyance(dry) = 0;
state = struct('stage_m', stage_m, 'conveyance_m3s', conveyance);
if ~face_properties
    return
end
area = horner(set.area_of_stage, set, index, interval, local);
hydrostatic = horner(set.hydrostatic_of_stage, set, index, interval, local);
next = breakpoint + set.count;
t = local ./ (set.stage_m(next) - set.stage_m(breakpoint));
top_width = (1 - t) .* set.top_width_m(breakpoint) + t .* set.top_width_m(next);
area(dry) = 0;
hydrostatic(dry) = 0;
top_width(dry) = 0;
state.area_m2 = area;
state.top_width_m = top_width;
state.hydrostatic_integral_m3 = hydrostatic;
end

function interval = find_interval(breakpoints, rows, index, query)
% Interval j with breakpoints(j) <= query < breakpoints(j+1), clamped to
% [1, rows-1], matching interp1/histc interval selection.
count = size(breakpoints, 1);
low = ones(numel(index), 1);
high = rows(index);
while true
    open = high - low > 1;
    if ~any(open)
        break
    end
    middle = floor((low + high) / 2);
    values = breakpoints(index + (middle - 1) * count);
    go_right = open & values(:) <= query;
    go_left = open & ~go_right;
    low(go_right) = middle(go_right);
    high(go_left) = middle(go_left);
end
interval = low;
end

function value = horner(coefficients, set, index, interval, local)
% Local pchip evaluation in the same Horner form as ppval and interp1.
plane = set.count * size(coefficients, 2);
base = index + (interval - 1) * set.count;
value = ((coefficients(base) .* local + coefficients(base + plane)) .* local + ...
    coefficients(base + 2 * plane)) .* local + coefficients(base + 3 * plane);
end

function [flux, maximum_speed] = hll_flux_vector( ...
        left, right, q_left, q_right, g, minimum_area)
count = numel(q_left);
u_left = zeros(count, 1);
u_right = zeros(count, 1);
c_left = zeros(count, 1);
c_right = zeros(count, 1);
wet_left = left.area_m2 > minimum_area;
wet_right = right.area_m2 > minimum_area;
u_left(wet_left) = q_left(wet_left) ./ left.area_m2(wet_left);
c_left(wet_left) = sqrt(g * left.area_m2(wet_left) ./ ...
    max(left.top_width_m(wet_left), eps));
u_right(wet_right) = q_right(wet_right) ./ right.area_m2(wet_right);
c_right(wet_right) = sqrt(g * right.area_m2(wet_right) ./ ...
    max(right.top_width_m(wet_right), eps));
physical_left = [q_left, q_left .* u_left + g * left.hydrostatic_integral_m3];
physical_right = [q_right, q_right .* u_right + g * right.hydrostatic_integral_m3];
left_speed = min(u_left - c_left, u_right - c_right);
right_speed = max(u_left + c_left, u_right + c_right);
maximum_speed = max(abs(left_speed), abs(right_speed));
flux = zeros(count, 2);
upwind_left = left_speed >= 0;
upwind_right = ~upwind_left & right_speed <= 0;
mixed = ~upwind_left & ~upwind_right & right_speed - left_speed > eps;
degenerate = ~upwind_left & ~upwind_right & ~mixed;
flux(upwind_left, :) = physical_left(upwind_left, :);
flux(upwind_right, :) = physical_right(upwind_right, :);
% Guarded so that scalar boundary calls never mix empty shapes.
if any(mixed)
    sl = left_speed(mixed);
    sr = right_speed(mixed);
    flux(mixed, :) = (sr .* physical_left(mixed, :) - ...
        sl .* physical_right(mixed, :) + ...
        sl .* sr .* ([right.area_m2(mixed), q_right(mixed)] - ...
        [left.area_m2(mixed), q_left(mixed)])) ./ (sr - sl);
end
if any(degenerate)
    flux(degenerate, :) = 0.5 * (physical_left(degenerate, :) + ...
        physical_right(degenerate, :));
end
end

% ---------------------------------------------------------------------------
% Scalar boundary solvers (same algorithms as HYDROHP_VARYING_STEP, with the
% table curves evaluated from precomputed coefficients)

function [state, stage_excess_m] = boundary_inner_state( ...
        face_section, cell_invert_m, cell_stage_m, boundary, side)
use_local_depth = side == "upstream" && ...
    isfield(boundary, 'upstream_state_type') && ...
    lower(string(boundary.upstream_state_type)) == "local_depth";
use_local_depth = use_local_depth || (side == "downstream" && ...
    isfield(boundary, 'downstream_type') && ...
    lower(string(boundary.downstream_type)) == "normal_depth");
if use_local_depth
    depth_m = max(cell_stage_m - cell_invert_m, 0);
    [state, stage_excess_m] = reconstructed_state_from_stage( ...
        face_section, face_section.invert_elevation_m + depth_m);
else
    [state, stage_excess_m] = reconstructed_state_from_stage( ...
        face_section, cell_stage_m);
end
end

function [state, stage_excess_m] = reconstructed_state_from_stage( ...
        section, stage_m)
stage_excess_m = max(stage_m - section.maximum_water_surface_elevation_m, 0);
if stage_excess_m > section.table_resolution_m
    state = state_from_stage(section, stage_m);
else
    state = state_from_stage(section, min( ...
        stage_m, section.maximum_water_surface_elevation_m));
end
end

function state = state_from_area(section, area_m2)
set = section.set;
maximum_area = set.maximum_area_m2;
area_tolerance = 1e-10 * max(1, maximum_area);
if area_m2 > maximum_area + area_tolerance
    error('hydrohp:AreaOutsideSectionRange', ...
        'Area %.6g m2 exceeds the local section range.', area_m2);
end
area_m2 = min(area_m2, maximum_area);
interval = find_interval(set.area_m2, set.rows, 1, area_m2);
depth_m = horner(set.depth_of_area, set, 1, interval, ...
    area_m2 - set.area_m2(interval));
state = state_from_stage(section, section.invert_elevation_m + depth_m);
end

function state = state_from_stage(section, stage_m)
state = stage_states(section.set, 1, stage_m, true);
end

function [flux, speed] = upstream_flux( ...
        section, inner, inner_discharge, prescribed_discharge, ...
        boundary, g, minimum_area)
if ~isscalar(prescribed_discharge) || ~isfinite(prescribed_discharge) || ...
        prescribed_discharge < 0
    error('hydrohp:InvalidUpstreamDischarge', ...
        'The prescribed upstream discharge must be a finite nonnegative scalar.');
end
if isfield(boundary, 'upstream_state_type') && ...
        lower(string(boundary.upstream_state_type)) == "local_depth"
    state = inner;
    if prescribed_discharge > 0
        critical = critical_state( ...
            section, prescribed_discharge, g, minimum_area);
        if state.area_m2 < critical.area_m2
            state = critical;
        end
    end
    flux = physical_flux(state, prescribed_discharge, g, minimum_area);
    speed = abs(safe_velocity( ...
        prescribed_discharge, state.area_m2, minimum_area)) + ...
        wave_celerity(state, g, minimum_area);
    return
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

function state = critical_state(section, discharge, g, minimum_area)
areas = section.table.area_m2;
widths = section.table.top_width_m;
wet = areas > minimum_area & widths > 0;
areas = areas(wet);
values = discharge^2 * widths(wet) - g * areas.^3;
upper_index = find(values <= 0, 1, 'first');
if isempty(upper_index)
    error('hydrohp:BoundaryStateUnavailable', ...
        'The critical inflow state exceeds the upstream section range.');
elseif upper_index == 1
    critical_area = areas(1);
else
    residual = @(area) discharge^2 * linear_curve( ...
        section.table.area_m2, section.table.top_width_m, area) - ...
        g * area.^3;
    critical_area = fzero(residual, ...
        [areas(upper_index - 1), areas(upper_index)]);
end
state = state_from_area(section, critical_area);
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
phi = zeros(size(area));
celerity = zeros(size(area));
positive = area > 0;
if any(positive)
    phi(positive) = pchip_curve(section.phi_of_area, area(positive));
    width = linear_curve(section.table.area_m2, section.table.top_width_m, ...
        area(positive));
    celerity(positive) = sqrt(g * area(positive) ./ max(width, eps));
end
end

function conveyance = conveyance_at_area(section, area)
conveyance = pchip_curve(section.conveyance_of_area, area);
end

function value = pchip_curve(curve, query)
% interp1(x, y, query, 'pchip'), including end-polynomial extrapolation.
breaks = curve.breaks;
count = numel(breaks);
interval = find_interval(breaks, count, ones(numel(query), 1), query(:));
local = query(:) - breaks(interval).';
c = curve.coefficients;
value = ((c(interval, 1) .* local + c(interval, 2)) .* local + ...
    c(interval, 3)) .* local + c(interval, 4);
value = reshape(value, size(query));
end

function value = linear_curve(x, y, query)
% interp1(x, y, query, 'linear'), NaN outside the table as in interp1.
count = numel(x);
interval = find_interval(x(:).', count, ones(numel(query), 1), query(:));
t = (query(:) - x(interval)) ./ (x(interval + 1) - x(interval));
value = (1 - t) .* y(interval) + t .* y(interval + 1);
value(query(:) < x(1) | query(:) > x(end)) = NaN;
value = reshape(value, size(query));
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
