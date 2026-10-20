function metrics = test_compiled_kernel_equivalence
%TEST_COMPILED_KERNEL_EQUIVALENCE Require identical compiled and reference steps.
%
% The compiled kernel is a performance rewrite only. Every state, flux
% diagnostic, and error it produces must equal the scalar reference kernel
% exactly, including wet/dry fronts, endpoint clamps, both boundary solvers,
% and complete runs with spin-up.

repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(repository_root);
compared_steps = 0;

% Lake at rest, clamped faces, and out-of-range errors.
[cells, faces, dx_m] = make_grid([10; 6; 10], [0; 0.4; 0], 40);
geometry = hydrohp_compile_geometry(cells, faces);
area = areas_at_stage(cells, 1.5 * ones(numel(cells), 1));
discharge = zeros(size(area));
compare_step(cells, faces, geometry, area, discharge, dx_m, 0.05, ...
    {'include_friction', false, 'cfl_target', 0.8});
clipped_faces = faces;
clipped_faces(2).maximum_water_surface_elevation_m = 1.5 - 0.001;
clipped = hydrohp_compile_geometry(cells, clipped_faces);
compare_step(cells, clipped_faces, clipped, area, discharge, dx_m, 0.05, ...
    {'include_friction', false, 'cfl_target', 0.8});
clipped_faces(2).maximum_water_surface_elevation_m = 1.5 - 0.03;
clipped = hydrohp_compile_geometry(cells, clipped_faces);
compare_error(cells, clipped_faces, clipped, area, discharge, dx_m, 0.05, ...
    {'include_friction', false}, 'hydrohp:StageOutsideSectionRange');
compare_error(cells, faces, geometry, area, discharge, dx_m, 100, ...
    {'include_friction', false, 'cfl_target', 0.7}, 'hydrohp:CFLViolation');
oversized = area;
oversized(5) = 2 * cells(5).table.area_m2(end);
compare_error(cells, faces, geometry, oversized, discharge, dx_m, 0.05, ...
    {}, 'hydrohp:AreaOutsideSectionRange');
compared_steps = compared_steps + 2;

% A wave through a contraction with friction, followed for 120 steps.
x = [cells.x_m].';
area = areas_at_stage(cells, 1.0 + 0.08 * exp(-((x - 25) / 8).^2));
discharge = zeros(size(area));
for step = 1:120
    [area, discharge] = compare_step(cells, faces, geometry, area, ...
        discharge, dx_m, 0.02, {'include_friction', true, 'cfl_target', 0.7});
    compared_steps = compared_steps + 1;
end

% A dam-break front into a dry bed exercises wet/dry cells and clipping.
[cells, faces, dx_m] = make_grid([8; 8; 8], [0; 0; 0], 60);
geometry = hydrohp_compile_geometry(cells, faces);
x = [cells.x_m].';
area = areas_at_stage(cells, 1.0 * (x < 40));
discharge = zeros(size(area));
adjustment_m3 = 0;
for step = 1:200
    [area, discharge, diagnostic] = compare_step(cells, faces, geometry, ...
        area, discharge, dx_m, 0.01, {'include_friction', true});
    adjustment_m3 = adjustment_m3 + diagnostic.volume_adjustment_m3;
    compared_steps = compared_steps + 1;
end

% Complete runs: normal-depth spin-up, hydrograph, steep inflow, and the
% characteristic upstream solver with a stage boundary.
sections = sloping_sections;
base = struct( ...
    'section_data', sections, ...
    'face_x_m', linspace(0, 100, 31).', ...
    'duration_s', 20, ...
    'downstream', struct('type', "normal_depth", 'slope', 0.001), ...
    'base_discharge_m3s', 1, ...
    'maximum_dt_s', 0.2, ...
    'spinup_duration_s', 100, ...
    'steady_tolerance_per_s', 1e-4, ...
    'hydrograph_time_s', [0; 10; 20], ...
    'hydrograph_discharge_m3s', [1; 1.5; 1]);
compare_run(base);
input = base;
input.section_data = steep_sloping_sections;
input.face_x_m = [0; 50; 100];
input.spinup_duration_s = 0;
compare_run(input);
input = base;
input.downstream = struct('type', "stage", 'stage_m', 0.8);
input.initialization_type = "stage";
input.initial_stage_m = 0.8;
input.base_discharge_m3s = 0;
input.spinup_duration_s = 0;
input.hydrograph_discharge_m3s = [0; 0.3; 0.1];
compare_run(input);

metrics = struct('compared_steps', compared_steps, ...
    'dry_bed_volume_adjustment_m3', adjustment_m3);
fprintf('HydroHP compiled-kernel equivalence checks passed (%d steps).\n', ...
    compared_steps);
end

function [area, discharge, diagnostic] = compare_step( ...
        cells, faces, geometry, area, discharge, dx_m, dt_s, options)
[area_ref, discharge_ref, diagnostic_ref] = hydrohp_varying_step( ...
    cells, faces, area, discharge, dx_m, dt_s, options{:});
[area, discharge, diagnostic] = hydrohp_varying_step_compiled( ...
    geometry, area, discharge, dx_m, dt_s, options{:});
assert(isequal(area, area_ref), 'Compiled area differs by %.3g m2.', ...
    max(abs(area - area_ref)));
assert(isequal(discharge, discharge_ref), ...
    'Compiled discharge differs by %.3g m3/s.', ...
    max(abs(discharge - discharge_ref)));
assert(isequaln(diagnostic, diagnostic_ref), 'Compiled diagnostics differ.');
end

function compare_error(cells, faces, geometry, area, discharge, dx_m, ...
        dt_s, options, identifier)
identifiers = strings(1, 2);
try
    hydrohp_varying_step(cells, faces, area, discharge, dx_m, dt_s, options{:});
catch exception
    identifiers(1) = exception.identifier;
end
try
    hydrohp_varying_step_compiled(geometry, area, discharge, dx_m, dt_s, ...
        options{:});
catch exception
    identifiers(2) = exception.identifier;
end
assert(all(identifiers == identifier), ...
    'Expected %s from both kernels, got %s and %s.', identifier, ...
    identifiers(1), identifiers(2));
end

function compare_run(input)
input.kernel = "reference";
reference = hydrohp_run_varying(input);
input.kernel = "compiled";
compiled = hydrohp_run_varying(input);
assert(compiled.metadata.kernel == "compiled");
reference.metadata = rmfield(reference.metadata, 'kernel');
compiled.metadata = rmfield(compiled.metadata, 'kernel');
% Face mass fluxes are reported by the compiled kernel only.
assert(isempty(reference.maxima.face_mass_flux_abs_m3s) && ...
    numel(compiled.maxima.face_mass_flux_abs_m3s) == numel(compiled.x_m) + 1);
reference.maxima = rmfield(reference.maxima, 'face_mass_flux_abs_m3s');
compiled.maxima = rmfield(compiled.maxima, 'face_mass_flux_abs_m3s');
reference.probes = rmfield(reference.probes, 'face_mass_flux_m3s');
compiled.probes = rmfield(compiled.probes, 'face_mass_flux_m3s');
assert(isequaln(compiled, reference), ...
    'Compiled and reference runs differ.');
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

function sections = sloping_sections
sections = trapezoid_sections([0.1; 0.05; 0]);
end

function sections = steep_sloping_sections
sections = trapezoid_sections([1; 0.5; 0]);
end

function sections = trapezoid_sections(invert_m)
chainage_m = [0; 50; 100];
rows = table();
for index = 1:numel(chainage_m)
    section = table( ...
        repmat("S" + index, 4, 1), repmat(chainage_m(index), 4, 1), ...
        [0; 0.8; 7.2; 8], invert_m(index) + [3; 0; 0; 3], ...
        0.03 * ones(4, 1), ...
        'VariableNames', {'section_id', 'chainage_m', 'offset_m', ...
        'elevation_m', 'manning_n'});
    rows = [rows; section]; %#ok<AGROW>
end
sections = hydrohp_build_section_tables(rows, 0.01, "scm");
end
