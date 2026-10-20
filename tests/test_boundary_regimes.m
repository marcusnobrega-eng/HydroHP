function metrics = test_boundary_regimes
%TEST_BOUNDARY_REGIMES P1: regime-based inflow and outlet boundaries.
%
% Inflow: the regime of the normal-depth state for the prescribed discharge
% decides whether one condition (subcritical: discharge, area from the
% outgoing characteristic) or two (supercritical: discharge and normal area)
% are imposed. Outlet: free outflow when supercritical, the normal-depth
% rating when a subcritical state exists, otherwise critical-depth control.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);

% Steep channel with a steep outlet: the legacy boundaries cannot close it.
steep = channel_input(trapezoid(8, 2, 4, 0.03), 0.03, 40);
legacy = steep;
legacy.boundary_formulation = "legacy";
identifier = "";
try
    hydrohp_run_varying(legacy);
catch exception
    identifier = string(exception.identifier);
end
assert(identifier == "hydrohp:BoundaryStateUnavailable");
result = hydrohp_run_varying(steep);
counts = result.initialization.ledger.boundary_steps;
steps = result.initialization.ledger.step_count;
assert(counts.upstream_supercritical_normal == steps);
assert(counts.downstream_rating == 0);
assert(counts.downstream_free_outflow + counts.downstream_critical_control == steps);
assert(result.initialization.ledger.wet_dry_adjustment_m3 == 0);
assert(result.initialization.relative_state_change_rate_per_s < 1e-6);
assert(result.metadata.boundary_formulation == "regime");

% Mild channel: one inflow condition and the normal-depth rating, with the
% same interior solution as the legacy boundaries.
mild = channel_input(trapezoid(8, 2, 4, 0.003), 0.003, 20);
result = hydrohp_run_varying(mild);
counts = result.initialization.ledger.boundary_steps;
steps = result.initialization.ledger.step_count;
assert(counts.upstream_subcritical_characteristic == steps);
assert(counts.downstream_rating == steps);
legacy = mild;
legacy.boundary_formulation = "legacy";
reference = hydrohp_run_varying(legacy);
inner = 5:numel(result.x_m) - 4;
relative_difference = max(abs(result.initialization.final_area_m2(inner) - ...
    reference.initialization.final_area_m2(inner)) ./ ...
    reference.initialization.final_area_m2(inner));
assert(relative_difference < 0.02, ...
    'Mild-channel interior differs from legacy by %.3g.', relative_difference);

% Mild channel ending on a steep outlet: critical-depth control.
outlet = mild;
outlet.downstream.slope = 0.05;
legacy = outlet;
legacy.boundary_formulation = "legacy";
identifier = "";
try
    hydrohp_run_varying(legacy);
catch exception
    identifier = string(exception.identifier);
end
assert(identifier == "hydrohp:BoundaryStateUnavailable");
result = hydrohp_run_varying(outlet);
counts = result.initialization.ledger.boundary_steps;
assert(counts.downstream_critical_control == result.initialization.ledger.step_count);
cells = hydrohp_map_sections_to_grid(outlet.section_data, result.x_m);
outlet_depth = hydrohp_depth_from_area(cells(end), ...
    result.initialization.final_area_m2(end));
[normal_depth, critical_depth] = reference_depths(cells(end), 2, 0.003);
assert(outlet_depth > critical_depth && outlet_depth < normal_depth, ...
    'Outlet depth %.3f m is not between critical %.3f and normal %.3f m.', ...
    outlet_depth, critical_depth, normal_depth);

% A hydrograph on an intermediate slope switches the inflow regime and uses
% every outlet control, without creating volume.
event = channel_input(trapezoid(8, 6, 6, 0.015), 0.015, 40);
event.hydrograph_time_s = [0; 300; 600; 1200];
event.hydrograph_discharge_m3s = [0.5; 60; 60; 0.5];
event.base_discharge_m3s = 0.5;
event.duration_s = 1200;
event.spinup_duration_s = 600;
event.maximum_dt_s = 2;
event.output_interval_s = 60;
result = hydrohp_run_varying(event);
counts = result.ledger.boundary_steps;
assert(counts.upstream_supercritical_normal > 0 && ...
    counts.upstream_subcritical_characteristic > 0);
assert(counts.downstream_free_outflow > 0 && counts.downstream_rating > 0 && ...
    counts.downstream_critical_control > 0);
assert(result.ledger.wet_dry_adjustment_m3 == 0);
assert(abs(result.ledger.relative_residual) < 1e-12);

% A sharp rise on a mild inlet admits no subcritical inflow state against
% the outgoing characteristic: the inflow chokes and enters at critical
% depth instead of stopping the run.
choke = channel_input(trapezoid(8, 6, 8, 0.003), 0.003, 40);
choke.hydrograph_time_s = [0; 10; 1200];
choke.hydrograph_discharge_m3s = [0.5; 80; 80];
choke.base_discharge_m3s = 0.5;
choke.duration_s = 1200;
choke.spinup_duration_s = 600;
choke.maximum_dt_s = 2;
choke.output_interval_s = 60;
choked = hydrohp_run_varying(choke);
assert(choked.status == "completed");
assert(choked.ledger.boundary_steps.upstream_two_condition_fallback > 0);
assert(choked.ledger.wet_dry_adjustment_m3 == 0);
assert(abs(choked.ledger.relative_residual) < 1e-12);

% A near-dry, motionless outlet cell (shallower than the first wet table
% level) leaves freely instead of stopping the step.
cells = hydrohp_map_sections_to_grid(mild.section_data, ...
    0.5 * (mild.face_x_m(1:end-1) + mild.face_x_m(2:end)));
faces = hydrohp_map_sections_to_grid(mild.section_data, mild.face_x_m);
geometry = hydrohp_compile_geometry(cells, faces);
area = zeros(numel(cells), 1);
for index = 1:numel(cells)
    area(index) = interp1(cells(index).table.depth_m, ...
        cells(index).table.area_m2, 0.3, 'pchip');
end
discharge = 2 * ones(size(area));
area(end) = interp1(cells(end).table.depth_m, cells(end).table.area_m2, 0.005, 'pchip');
discharge(end) = 0;
boundary = struct('upstream_discharge_m3s', 2, 'upstream_state_type', "local_depth", ...
    'downstream_type', "normal_depth", 'downstream_slope', 0.003);
[~, ~, diagnostic] = hydrohp_varying_step_compiled(geometry, area, discharge, ...
    mild.face_x_m(2) - mild.face_x_m(1), 0.1, boundary=boundary, ...
    boundary_formulation="regime", reconstruction="hydrostatic", ...
    return_face_mass_flux=true);
assert(diagnostic.downstream_boundary_code == 1);
assert(abs(diagnostic.face_mass_flux_m3s(end)) < 1e-3);

% Supercritical flow back toward the interior at the outlet reflects; the
% outlet never imports water.
area(end) = area(end - 1);
backflow = 2 * ones(size(area));
backflow(end) = -3 * area(end);    % u = -3 m/s, faster than c at 0.3 m depth
[~, ~, diagnostic] = hydrohp_varying_step_compiled(geometry, area, backflow, ...
    mild.face_x_m(2) - mild.face_x_m(1), 0.01, boundary=boundary, ...
    boundary_formulation="regime", reconstruction="hydrostatic", ...
    return_face_mass_flux=true);
assert(diagnostic.downstream_boundary_code == 4);
assert(abs(diagnostic.face_mass_flux_m3s(end)) < 1e-12);

% The reference kernel has only the legacy boundaries.
event.kernel = "reference";
event.reconstruction = "stage";
try
    hydrohp_run_varying(event);
    error('hydrohp:TestFailed', 'Regime boundaries ran on the reference kernel.');
catch exception
    assert(strcmp(exception.identifier, 'hydrohp:InvalidVaryingInput'));
end

metrics = struct('mild_relative_difference_from_legacy', relative_difference, ...
    'outlet_depth_m', outlet_depth, 'outlet_critical_depth_m', critical_depth, ...
    'outlet_normal_depth_m', normal_depth, 'event_boundary_steps', counts);
fprintf(['HydroHP boundary-regime checks passed: steep and steep-outlet cases ' ...
    'now close; event inflow %d supercritical / %d subcritical steps.\n'], ...
    counts.upstream_supercritical_normal, counts.upstream_subcritical_characteristic);
end

function input = channel_input(sections, outlet_slope, cell_count)
input = struct( ...
    'section_data', sections, ...
    'face_x_m', linspace(0, 600, cell_count + 1).', ...
    'hydrograph_time_s', [0; 1e5], ...
    'hydrograph_discharge_m3s', [2; 2], ...
    'duration_s', 1, ...
    'downstream', struct('type', "normal_depth", 'slope', outlet_slope), ...
    'base_discharge_m3s', 2, ...
    'spinup_duration_s', 1800, ...
    'steady_tolerance_per_s', 1e9, ...
    'maximum_dt_s', 5);
end

function [normal_depth, critical_depth] = reference_depths(section, discharge, slope)
data = section.table;
normal_area = interp1(data.conveyance_m3s, data.area_m2, ...
    discharge / sqrt(slope), 'pchip');
critical_area = fzero(@(area) discharge^2 * interp1(data.area_m2, ...
    data.top_width_m, area) - 9.81 * area.^3, [data.area_m2(2), data.area_m2(end)]);
normal_depth = hydrohp_depth_from_area(section, normal_area);
critical_depth = hydrohp_depth_from_area(section, critical_area);
end

function sections = trapezoid(bottom_width_m, side_m, bank_height_m, slope)
chainage_m = [0; 300; 600];
rows = table();
for index = 1:numel(chainage_m)
    invert = 20 - slope * chainage_m(index);
    section = table( ...
        repmat("S" + index, 4, 1), repmat(chainage_m(index), 4, 1), ...
        [0; side_m; side_m + bottom_width_m; 2 * side_m + bottom_width_m], ...
        invert + bank_height_m * [1; 0; 0; 1], 0.035 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section]; %#ok<AGROW>
end
sections = hydrohp_build_section_tables(rows, 0.02, "scm");
end
