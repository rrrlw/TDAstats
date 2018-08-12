context("Accurate visualization of known homologies")
library(TDAstats)

test_that("plot_persist produces ggplot object with correct data", {
  # use known circle2d dataset (provided with TDAstats)
  data("circle2d")
  circ.phom <- calculate_homology(circle2d)
  gpersist <- plot_persist(circ.phom)
  
  # make sure data in gpersist is correct for circle2d (expected values checked manually)
  #   extra as.character within as.integer to preserve factor levels as 0 and 1 (otherwise become 1 and 2)
  expect_equal(length(gpersist$data$dimension), 100)
  expect_equal(sum(as.integer(as.character(gpersist$data$dimension))), 1)
  expect_equal(sum(gpersist$data$birth == 0), 99)
  expect_equal(sum(gpersist$data$death > 0.21), 1)
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