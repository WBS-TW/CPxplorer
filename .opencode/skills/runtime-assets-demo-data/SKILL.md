---
name: runtime-assets-demo-data
description: Use when changing inst runtime assets, app instruction markdown/html, demo Skyline files, packaged Excel files, videos, images, or system.file paths.
---

# Runtime Assets Demo Data

Use this skill for packaged assets under `inst/`, including Shiny instruction files, demo data, Excel templates, videos, images, and paths loaded through `system.file(...)`.

## Repo Triggers

- `inst/instructions_CPions.md`
- `inst/instructions_CPions.html`
- `inst/instructions_CPquant.md`
- `inst/CPions_TP_formula.xlsx`
- `inst/cpxplorer-demo_2026-04-01/`
- `inst/CPxplorer_Logo.png`
- `system.file(...)`
- packaged demo or tutorial asset changes

## Runtime Rules

- Files under `inst/` are installed at the package root and loaded at runtime with `system.file(...)`.
- Preserve asset filenames and relative paths unless code references are updated at the same time.
- Demo Skyline Excel and `.sky*` files are integration fixtures for CPions and CPquant workflows.
- Avoid bloating package assets unless the new file is required for runtime, demos, or documentation.

## Verification

Choose checks based on the changed asset:

- `devtools::load_all(); CPions()` when CPions instructions/assets change
- `devtools::load_all(); CPquant()` when CPquant instructions/assets change
- `devtools::test(filter = "cpquant-demo-data")` when demo Skyline reports change
- `devtools::test(filter = "cpions")` when CPions demo/export fixtures change
- `devtools::build()` when package inclusion or installed paths matter

## Completion Checklist

- Runtime `system.file(...)` paths still resolve
- Relevant app instructions still render
- Demo data changes are covered by focused tests where possible
- Package build inclusion considered for new or renamed assets
