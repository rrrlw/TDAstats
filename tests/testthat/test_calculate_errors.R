context("Providing invalid parameters to calculate_homology")
library(TDAstats)

test_that("calculate_homology recognizes invalid parameters", {
  # mat is too small (not enough rows)
  test.mat <- matrix(1, nrow = 1, ncol = 3)
  expect_error(calculate_homology(test.mat),
               "Point cloud must have at least 2 points and at least 2 dimensions\\.")

  # mat is too small (not enough columns)
  test.mat <- matrix(1, nrow = 3, ncol = 1)
  expect_error(calculate_homology(test.mat),
               "Point cloud must have at least 2 points and at least 2 dimensions\\.")

  # mat has character elements
  test.mat <- matrix(c(1, 2, 3, "A"), nrow = 2, ncol = 2)
  expect_error(calculate_homology(test.mat),
               "Point cloud must contain values of class `numeric` or `integer` only\\.")

  # mat has missing elements
  test.mat <- matrix(c(1, 2, 3, NA), nrow = 2, ncol = 2)
  expect_error(calculate_homology(test.mat),
               "Point cloud has missing values\\.")

  # dim is numeric (not integer)
  test.mat <- matrix(1:4, nrow = 2, ncol = 2)
  expect_error(calculate_homology(test.mat, dim = 1.5),
               "dim parameter needs to be an integer")

  # dim is negative
  expect_error(calculate_homology(test.mat, dim = -1),
               "dim cannot be negative")

  # threshold is of type character
  expect_error(calculate_homology(test.mat, threshold = "hello"),
               "threshold parameter must be of type numeric")

  # invalid format parameter
  expect_error(calculate_homology(test.mat, format = "hello"),
               "format parameter should be either \"cloud\" or \"distmat\"")
})