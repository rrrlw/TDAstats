## Test environments

* local macOS: Sequoia (15.6.1) & Tahoe (26.2), R 4.5.2
  - `devtools::check()`
  - `devtools::check(env_vars = c('_R_CHECK_DEPENDS_ONLY_' = "true"))`
  - `devtools::check(manual = TRUE, remote = TRUE)`
* `win-builder`: old release, release, development

## R CMD check results

This resubmission follows a change of maintainer from Raoul Wadhwa to Cory Brunson.

Otherwise, there were no ERRORs, WARNINGs, or NOTEs.

