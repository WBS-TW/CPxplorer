---
name: cpions-skyline-interference
description: Use when changing CPions ion generation, interference detection, Skyline transition export, ion selection strategies, or related tests/snapshots.
---

# CPions Skyline Interference

Use this skill for changes involving CPions formula generation, isotope/adduct output, interference calculations, and Skyline transition-list export.

## Repo Triggers

- `R/CPions.R`
- `R/CPions_utils.R`
- `build_skyline_transition_list()`
- `compute_interference()`
- `compute_transition_interference()`
- `has_ms_interference()`
- Skyline Quan/Qual ion selection
- `tests/testthat/test-cpions-*.R`
- `tests/testthat/_snaps/cpions-*.md`

## Data Flow

- Normal and advanced ion tables are produced before interference analysis.
- Interference analysis populates state consumed later by the Skyline tab.
- Skyline export depends on the selected mode, quant ion setting, resolution, strategy, and preferred Qual-ion count.
- `CP_allions_compl2` is populated by the interference step and is expected by later CPions flows.

## Working Rules

- Keep CPions changes local to the affected generation, interference, or Skyline selection path.
- Preserve existing output column names unless the task explicitly changes the export contract.
- Treat `Note`, `Label Type`, `Ion Selection Rank`, `Selected Reason`, and `Interference at MS Res?` as user-visible export fields.
- Avoid editing the embedded `isotopes` table in `R/CPions_utils.R` unless the mass table itself is the intended change.
- Snapshot changes should be intentional and reviewed against the expected CPions output shape.

## Verification

Prefer focused testthat runs that match the changed path:

- `devtools::test(filter = "cpions-normal")`
- `devtools::test(filter = "cpions-advanced")`
- `devtools::test(filter = "cpions-interference")`
- `devtools::test(filter = "cpions-skyline")`
- `devtools::test(filter = "cpions-app-smoke")`

For app-level changes, also consider:

- `devtools::load_all(); CPions()`

## Completion Checklist

- Reviewed the upstream tab state needed by the changed CPions path
- Ran or attempted the relevant CPions focused tests
- Called out any snapshot changes explicitly
- Verified export columns when Skyline output changed
