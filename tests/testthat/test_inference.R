context("Test statistical inference on persistent homology")
library("TDAstats")

test_that("Basic permutation test (2-d) works correctly", {
  # set variables for reproducibility and maintenance ease
  set.seed(1)
  N.ITER <- 25
  
  # do permutation test on datasets provided with TDAstats
  data("circle2d")
  data("unif2d")
  perm.test <- permutation_test(circle2d, unif2d, iterations = N.ITER, update = 0)
  expect_equal(length(perm.test), 2)
  expect_equal(perm.test[[1]]$dimension, 0)
  expect_equal(perm.test[[2]]$dimension, 1)
  expect_equal(length(perm.test[[1]]$permvals), N.ITER)
  expect_equal(length(perm.test[[2]]$permvals), N.ITER)
  expect_true(perm.test[[1]]$pvalue < 0.05)
})

test_that("Basic permutation test (3-d) works correctly", {
  # set variables for reproducibility and maintenance ease
  set.seed(1)
  N.ITER <- 25
  
  # skip on CRAN b/c this will take longer
  skip_on_cran()
  data("sphere3d")
  data("unif3d")
  perm.test <- permutation_test(sphere3d, unif3d, iterations = N.ITER, update = 1, dim = 2)
  expect_equal(length(perm.test), 3)
  expect_equal(perm.test[[1]]$dimension, 0)
  expect_equal(perm.test[[2]]$dimension, 1)
  expect_equal(perm.test[[3]]$dimension, 2)
  expect_equal(length(perm.test[[1]]$permvals), N.ITER)
  expect_equal(length(perm.test[[2]]$permvals), N.ITER)
  expect_equal(length(perm.test[[3]]$permvals), N.ITER)
  expect_true(perm.test[[1]]$pvalue < 0.05)
})