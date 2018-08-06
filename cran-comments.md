## Test environments

* local Ubuntu 18.04 install, R 3.4.4 (old-release)
* Travis CI, R 3.5.0 (r-release)
* win-builder, R 3.5.0 (r-release)
* win-builder (r-devel)
* r-hub macOS 10.11, R 3.5.0 (r-release)

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE: the package was archived prior to the Aug 4 deadline, and appears as a new submission. The installation error has been fixed.

## Purpose

The primary purpose of such a quick update submission to CRAN is because of a clash with `index_t` in C++ code brought up by Dr. Brian Ripley that causes TDAstats to fail on Solaris systems (with an Aug 4 deadline to fix). This issue has been fixed in this update.
