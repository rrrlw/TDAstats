context("Accurate visualization of known homologies")
library(TDAstats)

test_that("plot methods detect invalid feature matrices", {
  # feature matrix with too few columns
  feature.mat <- matrix(1:6, ncol = 2)
  expect_error(plot_persist(feature.mat),
               "Invalid feature matrix with 2 columns: must have exactly 3 columns\\.")
  
  # feature matrix with too many columns
  feature.mat <- matrix(1:8, ncol = 4)
  expect_error(plot_barcode(feature.mat),
               "Invalid feature matrix with 4 columns: must have exactly 3 columns\\.")
  
  # feature matrix without any rows
  feature.mat <- matrix(NA, nrow = 0, ncol = 3)
  expect_error(plot_persist(feature.mat),
               "Invalid feature matrix: must have at least one row\\.")
  
  # feature matrix with first column containing characters
  feature.mat <- matrix(c("a", 2, 3, "b", 3, 4), ncol = 3, byrow = TRUE)
  expect_error(plot_barcode(feature.mat),
               "Invalid feature matrix: first column must be of type integer to represent feature dimension\\.")
  
  # feature matrix with second column containing characters (same result as third column containing characters)
  feature.mat <- matrix(c(0, "a", 2, 1, 2, 3), ncol = 3, byrow = TRUE)
  expect_error(plot_persist(feature.mat),
               "Invalid feature matrix: second and third columns must be of type numeric to represent feature birth and death\\.")

  # feature matrix with missing values
  feature.mat <- matrix(c(0, 1, 2, 1, 2, 3, 1, NA, 4), ncol = 3, byrow = TRUE)
  expect_error(plot_barcode(feature.mat),
               "Invalid feature matrix with 1 NAs: should have zero NAs.")
  
  # feature matrix with birth before death
  feature.mat <- matrix(c(0, 0, 1, 0, 0, 2, 1, 2, 1), ncol = 3, byrow = TRUE)
  expect_error(plot_persist(feature.mat),
               "Invalid feature matrix: birth must come before death for all features\\.")
})

test_that("plot_persist produces ggplot object with correct data", {
  # use known circle2d dataset (provided with TDAstats)
  data("circle2d")
  circ.phom <- calculate_homology(circle2d)
  gpersist <- plot_persist(circ.phom)
  flatpersist <- plot_persist(circ.phom, flat = TRUE)
  
  # make sure data in gpersist is correct for circle2d (expected values checked manually)
  #   extra as.character within as.integer to preserve factor levels as 0 and 1 (otherwise become 1 and 2)
  expect_equal(length(gpersist$data$dimension), 100)
  expect_equal(sum(as.integer(as.character(gpersist$data$dimension))), 1)
  expect_equal(sum(gpersist$data$birth == 0), 99)
  expect_equal(sum(gpersist$data$death > 0.21), 1)
  expect_true("ggplot" %in% class(flatpersist))
})

test_that("plot_barcode produces ggplot object with correct data", {
  # use known unif2d dataset (provided with TDAstats)
  data("unif2d")
  unif.phom <- calculate_homology(unif2d)
  gbarcode <- plot_barcode(unif.phom)
  
  # make sure data in gbarcode is correct for unif2d (expected values checked manually)
  #   extra as.character within as.integer to preserve factor levels as 0 and 1 (otherwise become 1 and 2)
  expect_equal(length(gbarcode$data$dimension), 117)
  expect_equal(sum(as.integer(as.character(gbarcode$data$dimension))), 18)
  expect_equal(sum(gbarcode$data$birth == 0), 99)
  expect_equal(sum(gbarcode$data$death > 0.25), 0)
})