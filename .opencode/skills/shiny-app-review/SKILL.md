---
name: shiny-app-review
description: Use when changing or reviewing the Shiny apps in R/CPions.R or R/CPquant.R. Focuses on reactive flow, tab dependencies, runtime assets, and app-specific verification.
---

# Shiny App Review

Use this skill for changes inside the Shiny app entrypoints or for bug reports involving app behavior, tabs, reactivity, uploads, or exports.

## Repo Triggers

Common triggers in this repo:
- `R/CPions.R`
- `R/CPquant.R`
- reactive bugs
- tab flow issues
- upload or export regressions
- app launch verification

## Repo Facts

- Most app logic lives directly inside the large `CPions()` and `CPquant()` functions.
- Shared helpers are in `R/CPions_utils.R` and `R/CPquant_utils.R`.
- App instructions are loaded from `inst/instructions_CPions.md` and `inst/instructions_CPquant.md`.

## CPions Constraints

- Interference and Skyline export flows depend on state built earlier in the app.
- `CP_allions_compl2` is populated by the interference step and is consumed later.
- When changing tab behavior, verify that earlier computed state is still available to later tabs.

## CPquant Constraints

- Treat `Group Mixtures` as the currently supported standards mode unless the task explicitly expands the model.
- Import logic expects Skyline Excel input and normalizes several possible Skyline column names.
- Changes to import logic should preserve existing Skyline column aliases.

## Review Priorities

- Reactive dependencies still fire in the right order.
- Objects created in one step remain available where later tabs expect them.
- Input validation still matches the current data model.
- Runtime file loads via `system.file(...)` still resolve to files under `inst/`.
- Download handlers and exports still receive the same data shape.

## Verification

Prefer targeted runtime verification:
- `devtools::load_all(); CPions()`
- `devtools::load_all(); CPquant()`

When possible, validate against demo files under `inst/cpxplorer-demo_2026-04-01/`.

## Change Style

- Keep edits local to the affected reactive path when possible.
- Avoid broad Shiny refactors unless the bug cannot be fixed locally.
- Add short comments only around non-obvious reactive or data-flow logic.

## Completion Checklist

- Confirm affected tab flow or import/export flow was reviewed end-to-end
- Verify app launch path relevant to the change
- Call out any runtime checks that could not be performed
