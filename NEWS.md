# v0.4.0

## Additions

* code coverage has been extended
* added support for distance matrix format for calculate_homology
* added `standardize` parameter and functionality to standardize point cloud size
* added `limit.num` parameter and functionality to consider only "significant" topological features

# v0.3.0

## Bug fixes

* fixed bug in `plot_persist` where setting a shorter axis limit would cause the reference diagonal to disappear (thank you, @corybrunson)

## Changes

* changing size of `plot_persist` output figures still preserves one-to-one length ratio of horizontal and vertical axes (fixed coordinate system; thank you, @corybrunson)

# v0.2.0

## Bug fixes

* fixed bug in `permutation_test` function that resulted in NAs in output
* fixed naming clash with `index_t` in C++ code; should make TDAstats compatible with Solaris systems

## Changes

* added `dim` parameter to `permutation_test` function to allow users to select maximum dimension for homology comparison

## Additions

* added automated testing using `testthat` package
* 4 sample datasets for users to learn/test TDAstats (and probably to follow future vignettes): unif2d, unif3d, circle2d, sphere3d
* added 2 vignettes introducing features of TDAstats to user through text and a case study
