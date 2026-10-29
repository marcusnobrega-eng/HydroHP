# HydroHP-1D

HydroHP is an open MATLAB solver of the unsteady Saint-Venant equations for
one-dimensional open-channel flow. It has two run paths.

| Path | Entry point | Geometry | Scheme |
|---|---|---|---|
| Varying-geometry kernel | `hydrohp_run_varying` | Surveyed cross-sections that vary along the reach | Finite volume, HLL fluxes, hydrostatic reconstruction |
| Legacy V5 | `hydrohp_run`, `HydroHP_1D_V5.m` | Regular sections, or one irregular section for the whole reach | Explicit Lax-Friedrichs |

The V5 path is the model published in Gomes Júnior et al. (2023). The
varying-geometry kernel extends it for cross-sections that vary along the
reach. Commit `9f1d944` is the kernel version used as the one-dimensional
baseline of a dam-break mapping study that is under review.

## Varying-geometry kernel

### What it adds

- Each surveyed cross-section gets its own hydraulic table. Area, wetted
  perimeter, top width, hydrostatic integral, and conveyance are tabulated by
  depth.
- The tables are mapped to every cell centre and cell face. Between two
  surveyed sections, each property is interpolated at the same fraction of the
  two sections' depth ranges. The bed is interpolated linearly.
- Each section is closed at the lower of its two bank tops. HydroHP never adds
  vertical walls. A flow that rises above a section top is reported, not
  extended.
- Every run keeps a mass ledger (inflow, outflow, storage, residual).

### Equations

Along chainage $x$, with no lateral inflow,

$$
\frac{\partial A}{\partial t}+\frac{\partial Q}{\partial x}=0,
$$

$$
\frac{\partial Q}{\partial t}
+\frac{\partial}{\partial x}\left(\frac{Q^2}{A}\right)
+gA\frac{\partial \eta}{\partial x}=-gAS_f,
\qquad
S_f=\frac{Q|Q|}{K^2},
$$

where $A$ is wetted area, $Q$ is discharge, $\eta$ is water-surface elevation,
$g$ is gravity, $S_f$ is friction slope, and $K$ is the section conveyance from
its table. For one channel with Manning coefficient $n$ and hydraulic radius
$R$, $K = AR^{2/3}/n$, so $S_f = n^2Q|Q|/(A^2R^{4/3})$.

### Numerics

- Finite-volume scheme on a uniform grid with HLL fluxes (Harten et al., 1983)
  at every face. Friction is applied semi-implicitly after the flux update.
- Hydrostatic reconstruction (Audusse et al., 2004) is the runner default. At
  each interior face the bed is the higher of the two neighbouring cell
  inverts, and each side's depth is `max(0, eta - face bed)`. The
  reconstructed depth never exceeds the cell depth, so depths stay nonnegative
  on dry beds and at wet and dry fronts. Equal water levels give equal face
  states, so a lake at rest stays at rest.
- Stage reconstruction (`reconstruction = "stage"`) is kept only for the tests
  that require the compiled kernel to match the reference kernel exactly.
- The time step follows a CFL target (default 0.7) between `minimum_dt_s` and
  `maximum_dt_s`. A step that violates the CFL condition is halved and
  repeated.
- Two kernels give identical results for the same options. `"compiled"`
  (default) is vectorized and is required for hydrostatic reconstruction and
  regime boundaries. `"reference"` is the scalar version kept for the
  equivalence tests.

### Boundaries

With `boundary_formulation = "regime"` (default), both boundaries follow the
local flow regime.

Upstream, the discharge hydrograph enters as follows.

| Inlet state | Rule |
|---|---|
| Normal-depth state of the inflow is supercritical | Impose the discharge and its normal-depth area |
| Subcritical | Impose the discharge. Take the area from the outgoing characteristic |
| No outgoing characteristic (dry inlet cell, or supercritical cell under subcritical normal flow), or no subcritical state admits the inflow | Impose the discharge at normal depth if the inlet cell is dry, at critical depth otherwise |
| Zero discharge | Reflective wall |

Downstream, the outlet works as follows.

| Outlet state | Rule |
|---|---|
| Supercritical, dry, or nearly dry | Free outflow |
| Subcritical | Normal-depth rating at `downstream.slope`, or a fixed `downstream.stage_m` |
| Subcritical cell but no subcritical normal-depth state (steep outlet) | Critical-depth control |
| Supercritical flow back into the reach | Reflected, so the outlet never imports water |

The outlet rating uses one constant energy slope for the whole run. HydroHP
does not compute it. The caller supplies it as `downstream.slope`.

`boundary_formulation = "legacy"` keeps the earlier inlet, which extrapolates
the local depth and falls back to critical depth.

### Initial state

| `initialization_type` | Start |
|---|---|
| `"normal_depth"` (default) | Normal depth of `base_discharge_m3s` on the local bed slope. Pools and very flat cells use the reach-average slope. Add a spin-up with `spinup_duration_s` |
| `"dry"` | Zero area and zero discharge. The hydrograph alone wets the reach. No spin-up is needed |
| `"stage"` | A given water level, one value or one per cell |

The dam-break study used the dry start.

### Applicability guard

When water rises above the top of a section, the run cannot continue on that
geometry. With the default `overtopping_action = "error"`, HydroHP raises an
error. With `"stop"`, the run ends at the last accepted state,
`result.status` is `"overtopped"`, and `result.overtopping` records the time
and chainage. This is a limit of the section geometry, not a solver failure.

### Run it

Run from the repository root. This example uses the three sections in
`Examples/cross_sections_example.csv` and starts dry.

```matlab
addpath(pwd)
section_data = hydrohp_build_section_tables( ...
    "Examples/cross_sections_example.csv", 0.02, "scm");

input = struct( ...
    'section_data', section_data, ...
    'face_x_m', linspace(0, 1000, 41).', ...   % 40 cells of 25 m
    'hydrograph_time_s', [0; 600; 1800; 3600; 7200], ...
    'hydrograph_discharge_m3s', [0; 40; 15; 0; 0], ...
    'duration_s', 7200, ...
    'downstream', struct('type', "normal_depth", 'slope', 1e-3), ...
    'initialization_type', "dry", ...
    'output_interval_s', 60);
result = hydrohp_run_varying(input);

result.status                    % "completed" or "overtopped"
result.maxima.stage_m            % maximum water level per cell (m)
result.maxima.velocity_m_s       % maximum mean velocity per cell (m/s)
result.ledger.relative_residual  % volume balance
```

`hydrohp_postprocess_varying(result, input)` writes the V5 post-processing
figures, videos, and CSV files for the run, with the same names and format.
`Examples/run_varying_example.m` runs this case and writes them.
The [user manual](docs/user_manual.md) walks through a full run, from
sections to checks and troubleshooting.
[docs/inputs_outputs.md](docs/inputs_outputs.md) lists every input, option,
and output.

## Legacy V5 path

Run a workbook without editing the solver:

```matlab
input = hydrohp_read_excel("Examples/Example_2_Dam_Break.xlsx");
result = hydrohp_run(input);
```

The returned structure contains time, distance, area, discharge, depth,
velocity, Froude, Courant, and execution-time arrays. Set
`input.postprocess = true` to request the legacy figures and CSV outputs.
`Examples/run_from_script.m` also shows optional script-defined hydrograph and
Manning values. The original `HydroHP_1D_V5.m` entry point remains available.

The example workbooks cover regular channel design, a regular-section dam
break, and a single irregular section. The V5 path supports discharge or stage
hydrographs, trapezoidal, circular, parabolic, and irregular sections, single
and divided-channel hydraulic properties, and automatic plots, videos, and CSV
outputs.

## Requirements

- MATLAB. The frozen V5 baseline was run with R2023b. The test suite passes
  with R2025b.
- No toolbox is needed for the varying-geometry kernel.
- Excel-compatible `.xlsx` files for the V5 inputs.

## Tests

Each file in `tests/` is a function that stops with an error when a check
fails. Run all of them from the repository root:

```bash
matlab -batch "addpath(pwd); addpath('tests'); files = dir('tests/test_*.m'); for k = 1:numel(files), [~, name] = fileparts(files(k).name); disp(name); feval(name); end"
```

| Test | Checks |
|---|---|
| `test_cross_section_input` | Section table from MATLAB, CSV, and Excel, ordering, and input errors |
| `test_section_hydraulics` | Section tables against analytic geometry |
| `test_section_grid_mapping` | Bracketing and interpolation of tables between sections |
| `test_varying_geometry_kernel` | Lake at rest, CFL control, wet and dry fronts, volume |
| `test_compiled_kernel_equivalence` | Compiled and reference kernels give identical steps |
| `test_hydrostatic_reconstruction` | Positivity and balance of hydrostatic reconstruction |
| `test_negative_area_capture` | Capture of the first negative raw area on a steep reach |
| `test_varying_run` | Runner boundaries, initialization, and mass ledger |
| `test_boundary_regimes` | Regime-based inflow and outlet rules |
| `test_pool_initialization` | Normal-depth start on a bed with pools |
| `test_dry_start` | Dry start, exact inflow volume, volume balance |
| `test_event_outputs` | Sparse history, online maxima, and probes |
| `test_overtopping_guard` | Stop and report when water leaves a section |
| `test_postprocess_varying` | Every post-processing file, CSV values, probe and history sources |
| `test_h7_regression_package` | Published normal flow and analytical dam break |
| `test_hydrohp_interface` | Legacy V5 interface from Excel and from a structure |

## Repository layout

| Path | Contents |
|---|---|
| `hydrohp_*.m` | Varying-geometry kernel, section tables, grid mapping, post-processing, and the callable V5 interface |
| Other root `.m` files | V5 solver, input, hydraulic-property, and postprocessing scripts |
| `Examples/` | Example inputs for both paths |
| `tests/` | Test functions |
| `docs/user_manual.md` | Step-by-step guide to the varying-geometry kernel |
| `docs/inputs_outputs.md` | Inputs, options, and outputs of the varying-geometry kernel |
| `docs/cross_section_input.md` | Cross-section table schema and validation rules |
| `docs/reference/` | Published supplement and hydraulic-property reference figures |
| `docs/tutorial/` | Tutorial slides and teaching material |
| `apps/` | Optional MATLAB applications |

Generated results and historical ZIP distributions are not kept on the active
branch. They remain recoverable from Git history and the tag
`legacy-v5-pre-modernization-2026-10-07`.

## Citation

Please cite the model paper:

> Gomes Júnior, M. N., Castillo Rápalo, L. M. C., Oliveira, P. T. S.,
> Giacomoni, M. H., do Lago, C. A. F., and Mendiondo, E. M. (2023).
> Modeling unsteady and steady 1D hydrodynamics under different hydraulic
> conceptualizations: Model/software development and case studies.
> *Environmental Modelling & Software*, 167, 105733.
> https://doi.org/10.1016/j.envsoft.2023.105733

If you use the varying-geometry kernel, also give the commit you ran. The
dam-break study used commit `9f1d944`. Machine-readable citation metadata are
in `CITATION.cff`.

Numerical methods:

- Harten, A., Lax, P. D., and van Leer, B. (1983). On upstream differencing
  and Godunov-type schemes for hyperbolic conservation laws. *SIAM Review*,
  25(1), 35-61. https://doi.org/10.1137/1025002
- Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., and Perthame, B.
  (2004). A fast and stable well-balanced scheme with hydrostatic
  reconstruction for shallow water flows. *SIAM Journal on Scientific
  Computing*, 25(6), 2050-2065. https://doi.org/10.1137/S1064827503431090

## Known limitations

- The varying-geometry runner needs a uniform grid that lies between the first
  and last section chainages.
- There is no lateral inflow.
- Velocity is the mean over each cross-section. Near a dry front its maximum
  can come from water only a few centimetres deep.
- Hydrostatic reconstruction overestimates depth where the bed drop across a
  face exceeds the water depth (Delestre et al., 2012). The kernel counts
  those faces in its step diagnostics.
- The bed between sections is interpolated linearly, so it can differ from the
  true bed on long section spacings.
- The V5 kernel remains a monolithic script and still reads its base
  configuration and geometry from an Excel workbook. Some legacy
  postprocessing profiles fail after the numerical loop completes.

## License

HydroHP is distributed under the [MIT License](LICENSE).

## Contact

Marcus Nóbrega Gomes Júnior, `marcusnobrega.engcivil@gmail.com`
