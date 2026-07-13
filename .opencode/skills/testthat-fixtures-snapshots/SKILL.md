---
name: testthat-fixtures-snapshots
description: Use when changing tests/testthat, fixture helpers, snapshot files, or test verification commands for this R package.
---

# Testthat Fixtures Snapshots

Use this skill for test changes, fixture updates, snapshot expectation changes, or verification planning involving this package's `tests/testthat` suite.

## Repo Triggers

- `tests/testthat.R`
- `tests/testthat/*.R`
- `tests/testthat/helper-*.R`
- `tests/testthat/_snaps/*.md`
- `expect_snapshot_value()`
- fixture helpers such as `normal_fixture()`, `skyline_fixture()`, and `cpquant_demo_base_df()`

## Test Suite Shape

- The package uses testthat edition 3.
- CPions has focused normal, advanced, interference, Skyline, and Shiny server smoke tests.
- CPquant has demo-data tests covering GC and LC Skyline reports, calibration, isotope QC, and deconvolution helpers.
- Snapshot files under `_snaps` are generated artifacts of test expectations, but they are still reviewed source-control artifacts.

## Working Rules

- Prefer focused `devtools::test(filter = "...")` runs while iterating.
- Run broader tests when shared helpers or package-wide behavior changes.
- Do not update snapshots just to make tests pass; first confirm the changed output is intentional.
- Keep fixture helpers small and deterministic.
- When snapshots change, summarize what changed in the serialized output shape or values.

## Useful Commands

- `devtools::test()`
- `devtools::test(filter = "cpions-skyline")`
- `devtools::test(filter = "cpquant-demo-data")`
- `testthat::snapshot_review()`
- `testthat::snapshot_accept()`

## Completion Checklist

- Relevant focused tests were run or attempted
- Snapshot changes were reviewed and explained
- Fixture changes preserve deterministic output
- Final response states any tests that could not be run
