# TDAstats 0.2.0

## Bug fixes

* fixed bug in `permutation_test` function that resulted in NAs in output
* fixed naming clash with `index_t` in C++ code; should make TDAstats compatible with Solaris systems

## Changes

* added `dim` parameter to `permutation_test` function to allow users to select maximum dimension for homology comparison

## Additions

* added automated testing using `testthat` package
