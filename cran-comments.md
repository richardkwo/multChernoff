## Resubmission

This is a new submission.

## Test environments

- Local macOS Sonoma 14.6.1, R 4.5.2
- Additional platforms: not yet run

## R CMD check results

- Current local tarball check under UTF-8 locale completed without code,
  documentation, or test failures.
- Remaining local issues were environment-specific:
  - vignette HTML was not built because `pandoc` is not available on this
    machine
  - CRAN incoming URL checks could not complete because network access was not
    available in the current environment
- Replace this section with the final clean `R CMD check --as-cran` results
  before submission.

## Comments

- The package provides finite-sample tail bounds and critical values for
  multinomial likelihood ratio tests.
- The submission includes package-level tests and a lightweight vignette.
- Before final submission, update the test-environment list above with all
  platforms used for the final checks, such as win-builder or R-hub.
