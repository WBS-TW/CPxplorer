---
name: cran-release-check
description: Use when changing DESCRIPTION, .Rbuildignore, cran-comments.md, CRAN-SUBMISSION, package versioning, R CMD check, or release prep.
---

# CRAN Release Check

Use this skill for package metadata, build configuration, CRAN notes, release preparation, and package-check workflows.

## Repo Triggers

- `DESCRIPTION`
- `.Rbuildignore`
- `cran-comments.md`
- `CRAN-SUBMISSION`
- package version changes
- dependency changes
- `R CMD check`
- `devtools::check()`
- release notes or CRAN submission requests

## Working Rules

- Keep `DESCRIPTION` fields valid for R packaging tools.
- Update generated documentation with `devtools::document()` when exports or roxygen comments change.
- Keep `.Rbuildignore` aligned with files that should not ship in the package tarball.
- Do not claim CRAN readiness without running or attempting package checks.
- Record check notes accurately in `cran-comments.md` only after verifying current results.

## Verification

Use the narrowest check that matches the task, escalating as needed:

- `devtools::document()` for roxygen/export changes
- `devtools::test()` for package test coverage
- `devtools::check()` or `rcmdcheck::rcmdcheck()` for release readiness

## Completion Checklist

- Metadata remains valid and internally consistent
- Generated docs refreshed if needed
- Tests/checks run or attempted and results reported
- CRAN notes match actual check output
