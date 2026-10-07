# HydroHP-1D

HydroHP-1D is a MATLAB research model for steady and unsteady one-dimensional
open-channel flow. It solves the Saint-Venant equations with an explicit
Lax-Friedrichs formulation and supports channel-design, flood-routing, and
dam-break applications.

## Current capabilities

- prescribed discharge or stage hydrographs;
- trapezoidal, circular, parabolic, and irregular cross-sections;
- single-channel and divided-channel hydraulic-property methods;
- spatially varying bed elevation, roughness, and regular-section dimensions;
- automatic plots, videos, and CSV outputs.

The modernization branch also includes a validated long-table input format
and hydraulic-property tables for spatially varying irregular cross-sections.
The numerical kernel does not use those profiles yet.

The current V5 irregular-section path uses one profile throughout the modeled
reach. Support for several spatially varying irregular profiles is under
development.

## Requirements

- MATLAB; the frozen V5 baseline was executed with R2023b;
- Excel-compatible `.xlsx` input files.

## Quick start

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

The example inputs cover regular channel design, a regular-section dam break,
and a single irregular section.

## Repository layout

| Path | Contents |
|---|---|
| Root `.m` files | V5 solver, callable interface, input, hydraulic-property, and postprocessing scripts |
| `Examples/` | Runnable Excel inputs |
| `tests/` | Callable-interface and input-schema checks |
| `docs/cross_section_input.md` | Multi-section table schema and validation rules |
| `docs/reference/` | Published supplement and hydraulic-property reference figures |
| `docs/tutorial/` | Tutorial slides and teaching material |
| `apps/` | Optional MATLAB applications |

Generated results and historical ZIP distributions are not kept on the active
branch. They remain recoverable from Git history and the tag
`legacy-v5-pre-modernization-2026-10-07`.

## Citation

Please cite:

> Gomes Júnior, M. N., Castillo Rápalo, L. M. C., Oliveira, P. T. S.,
> Giacomoni, M. H., do Lago, C. A. F., and Mendiondo, E. M. (2023).
> Modeling unsteady and steady 1D hydrodynamics under different hydraulic
> conceptualizations: Model/software development and case studies.
> *Environmental Modelling & Software*, 167, 105733.
> https://doi.org/10.1016/j.envsoft.2023.105733

Machine-readable citation metadata are provided in `CITATION.cff`.

## Known limitations

- The V5 numerical kernel remains a monolithic script and still obtains its
  base configuration and geometry from an Excel workbook.
- Several varying irregular cross-sections are not yet supported.
- Automated regression tests are being added during modernization.
- Some legacy postprocessing profiles fail after the numerical loop completes.

## License

HydroHP is distributed under the [MIT License](LICENSE).

## Contact

Marcus Nóbrega Gomes Júnior — `marcusnobrega.engcivil@gmail.com`
