# Varying-geometry kernel: inputs and outputs

This page describes `hydrohp_run_varying`. The legacy V5 path (`hydrohp_run`)
keeps its Excel workbook format.

## Inputs

A run takes one MATLAB structure. Six fields are required.

| Field | Meaning |
|---|---|
| `section_data` | Section tables from `hydrohp_build_section_tables` |
| `face_x_m` | Cell-face chainages (m). Uniform spacing, at least three faces, inside the surveyed chainage range |
| `hydrograph_time_s` | Hydrograph times (s). Starts at 0, strictly increasing, reaches `duration_s` |
| `hydrograph_discharge_m3s` | Inflow at those times (m3/s). Nonnegative |
| `duration_s` | Event length (s) |
| `downstream` | Outlet boundary, see below |

### Section geometry

Sections arrive as one long table with one row per profile point. The same
columns work in a MATLAB table, a CSV file, or the `CrossSections` sheet of an
Excel workbook.

| Column | Meaning |
|---|---|
| `section_id` | Profile identifier, shared by all its points |
| `chainage_m` | Distance along the reach (m). One value per profile |
| `offset_m` | Transverse coordinate (m) |
| `elevation_m` | Ground elevation (m), one vertical datum for all profiles |
| `manning_n` | Manning coefficient of the segment that starts at the point |
| `x_m`, `y_m` | Optional plan coordinates |

```text
section_id,chainage_m,offset_m,elevation_m,manning_n
XS001,0,0,105.0,0.060
XS001,0,20,100.0,0.035
XS001,0,40,105.0,0.060
XS002,500,0,104.2,0.060
...
```

[cross_section_input.md](cross_section_input.md) gives the full validation
rules. Build the tables with

```matlab
section_data = hydrohp_build_section_tables(source, depth_step_m, method);
```

where `source` is a table, CSV path, or workbook path, `depth_step_m` is the
largest depth increment of the table (surveyed vertex elevations are also
inserted), and `method` is `"scm"` (one channel with an equivalent roughness)
or `"dcm"` (conveyance summed over roughness zones). Each table stops at the
lower of the two outer-bank elevations. HydroHP does not add vertical walls,
so a section must be wide enough for the flow it will carry.

HydroHP takes the sections as given. It does not decide where they go or how
wide they are. The section rule used in the dam-break study (target stage
above the larger of the critical and normal stages of the peak inflow) lives
in the study's own repository, in `scripts/preflight_hydrohp_sections_v4.py`.

### Hydrograph

`hydrograph_time_s` and `hydrograph_discharge_m3s` are column vectors of the
same length. HydroHP interpolates them linearly and never steps across a
hydrograph breakpoint. A hydrograph longer than the event is fine. Only the
part up to `duration_s` is used.

### Downstream boundary

| `downstream` | Outlet |
|---|---|
| `struct('type', "normal_depth", 'slope', S)` | Normal-depth rating with one constant energy slope `S` |
| `struct('type', "stage", 'stage_m', Z)` | Fixed water level `Z` |

Under the default regime boundaries, a supercritical or dry outlet drains
freely whatever the type. The normal-depth start requires a `"normal_depth"`
outlet with a positive slope. With a `"stage"` outlet, start dry or from a
stage.

### Options

All options have defaults.

| Field | Default | Meaning |
|---|---|---|
| `initialization_type` | `"normal_depth"` | `"dry"`, `"stage"`, or `"normal_depth"` |
| `base_discharge_m3s` | first hydrograph value | Discharge of the normal-depth start and spin-up |
| `initial_stage_m` | `NaN` | Water level for the `"stage"` start, one value or one per cell |
| `spinup_duration_s` | `0` | Spin-up length at the base discharge before the event |
| `steady_initialization` | `"require"` | Stop if the spin-up is not steady, or `"report"` it and run the event |
| `steady_tolerance_per_s` | `1e-4` | Relative change rate that counts as steady |
| `reconstruction` | `"hydrostatic"` | `"stage"` only for kernel-equivalence tests |
| `boundary_formulation` | `"regime"` | `"legacy"` keeps the earlier inlet |
| `kernel` | `"compiled"` | `"reference"` is the scalar kernel for tests |
| `cfl_target` | `0.7` | CFL number of the adaptive time step |
| `minimum_dt_s`, `maximum_dt_s` | `1e-5`, `1` | Time-step bounds (s) |
| `minimum_area_m2` | `1e-10` | Area below which a cell is dry |
| `include_friction` | `true` | Friction source term |
| `gravity_m_s2` | `9.81` | Gravity |
| `overtopping_action` | `"error"` | `"stop"` ends the run when water leaves a section |
| `output_interval_s` | `0` | Spacing of the stored history (s). `0` stores every step |
| `probe_x_m` | empty | Chainages for per-step series |
| `progress_file` | `""` | Tab-separated progress log |
| `checkpoint_file` | `""` | `.mat` file with the latest state |
| `diagnostic_wall_interval_s` | `300` | Wall-clock spacing of progress and checkpoint writes (s) |

## Outputs

All values are per cell of the solver grid, at the cell centres
`result.x_m`, not at the surveyed sections. Velocity is the mean over the
cross-section.

| Field | Contents |
|---|---|
| `status` | `"completed"` or `"overtopped"` |
| `overtopping` | Phase, model time, step, chainage, and message of an overtopping stop |
| `x_m` | Cell-centre chainages (m) |
| `time_s` | Times of the stored history (s) |
| `area_m2`, `discharge_m3s` | History, one row per stored time, one column per cell |
| `final_area_m2`, `final_discharge_m3s` | State at the end of the run |
| `maxima` | Maxima over every time step, see below |
| `probes` | Per-step series at `probe_x_m` (area, discharge, depth, face flux) |
| `ledger` | Volume balance and step statistics |
| `initialization` | Start type, spin-up result, and slope fallbacks |
| `metadata` | Kernel, reconstruction, boundaries, and time-step settings used |

### Maxima

The maxima are updated at every accepted step, so a sparse history does not
change them.

| Field | Contents |
|---|---|
| `maxima.stage_m` | Maximum water level (m) |
| `maxima.depth_m` | Maximum depth above the cell invert (m) |
| `maxima.velocity_m_s` | Maximum mean velocity (m/s) |
| `maxima.depth_velocity_m2_s` | Maximum of depth times velocity (m2/s) |
| `maxima.area_m2`, `maxima.discharge_abs_m3s` | Maximum area and absolute discharge |
| `maxima.time_of_maximum_depth_s`, `maxima.time_of_maximum_velocity_s` | When the maxima occurred (s) |
| `maxima.face_mass_flux_abs_m3s` | Maximum absolute flux at each face |

### Stage, depth, and velocity histories

The history stores area and discharge. Convert them with the same mapped
section tables the solver used:

```matlab
cell_sections = hydrohp_map_sections_to_grid(input.section_data, result.x_m);
depth_m = zeros(size(result.area_m2));
for k = 1:numel(cell_sections)
    depth_m(:, k) = hydrohp_depth_from_area(cell_sections(k), result.area_m2(:, k));
end
stage_m = depth_m + [cell_sections.invert_elevation_m];
velocity_m_s = zeros(size(result.area_m2));
wet = result.area_m2 > 1e-10;
velocity_m_s(wet) = result.discharge_m3s(wet) ./ result.area_m2(wet);
```

### Post-processing

`files = hydrohp_postprocess_varying(result, input, Name, Value)` writes the
V5 figures, videos, and CSV files. The user manual lists the products. Every
option has a default.

| Option | Default | Meaning |
|---|---|---|
| `output_folder` | `"Modeling_Results"` | Folder for all files, created if missing |
| `run_id`, `run_name` | `"1"`, `"HydroHP"` | File name prefix `<run_id>_<run_name>_` |
| `observation_x_m` | probes, or inlet, middle, outlet | Up to five chainages |
| `observation_labels` | `"x = ... m"` | One legend label per point |
| `figures`, `videos`, `csv` | `true` | Which product groups to write |
| `csv_layout` | `"by_section"` | `"by_state"` groups `Detailed_Output.csv` by state |
| `arrival_depth_m` | `0.5` | Depth rise above the start that marks flood arrival |
| `person_height_m`, `person_width_m` | `1.70`, `0.35` | Person used for the V5 force |
| `specific_weight_n_m3` | `9810` | Specific weight of water |
| `computational_time_s` | `NaN` | Wall time to report in `Dam_Breach_Like_Outputs.csv` |
| `frame_rate` | `20` | Video frames per second |

### Ledger

| Field | Contents |
|---|---|
| `cumulative_inflow_m3`, `cumulative_outflow_m3` | Boundary volumes |
| `initial_storage_m3`, `final_storage_m3` | Reach storage |
| `wet_dry_adjustment_m3` | Volume changed by clipping negative areas to zero |
| `residual_m3`, `relative_residual` | Storage change minus net inflow |
| `step_count`, `minimum_accepted_dt_s`, `maximum_cfl_number` | Step statistics |
| `stage_endpoint_clamp_count`, `maximum_stage_endpoint_excess_m` | Face stages clamped at a table end |
| `boundary_steps` | Steps spent in each inflow and outlet rule (regime boundaries only) |
