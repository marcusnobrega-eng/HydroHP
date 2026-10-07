# Spatially varying cross-section input

HydroHP uses one long table for surveyed cross-sections. Each row is one
profile point. The same columns are used in MATLAB tables, CSV files, and the
`CrossSections` sheet of an Excel workbook.

| Column | Requirement |
|---|---|
| `section_id` | Nonempty identifier shared by all points in one profile |
| `chainage_m` | Nonnegative longitudinal position; one unique value per profile |
| `offset_m` | Transverse profile coordinate; no duplicates within a profile |
| `elevation_m` | Elevation in metres in one consistent vertical datum |
| `manning_n` | Positive Manning coefficient for the segment beginning at the point |
| `x_m`, `y_m` | Optional plan coordinates; provide both or neither |

Each profile needs at least three points. Rows may arrive in any order because
the identifiers and coordinates define their position. The reader returns a
canonical order by `chainage_m`, `section_id`, and `offset_m` and records a
warning when it had to reorder the input.

```matlab
[sections, report] = hydrohp_read_cross_sections( ...
    "Examples/cross_sections_example.csv");
```

The function also accepts a MATLAB table or an Excel workbook containing a
sheet named `CrossSections`. It stops with one combined error message when the
table is invalid. `report` records the number of profiles and points, sorted
chainages, errors, and warnings.

The H3 reader validates structure only. Hydraulic-property tables, valid stage
ranges, and interpolation between profiles are implemented and tested in H4
and H5.
