context("Make sure different input formats result in same output")
library("TDAstats")

# setup functions
# calculates the distance between two points
calc.dist <- function(point1, point2) {
  sqrt(sum((point1 - point2) ^ 2))
}

# calculates a distance matrix for a point cloud
calc.distmat <- function(point.cloud) {
  # create empty matrix
  ans.mat <- matrix(NA, nrow = nrow(point.cloud), ncol = nrow(point.cloud))
  
  # populate matrix
  for (i in 1:nrow(point.cloud)) {
    for (j in 1:nrow(point.cloud)) {
      ans.mat[i, j] <- calc.dist(point.cloud[i, ], point.cloud[j, ])
    }
  }
  
  # return distance matrix
  return(ans.mat)
}

# only checks 2-d data (3-d too long, maybe skip_cran() on that in a separate test)
test_that("Point cloud and lower distance matrix formats are equivalent", {
  # generate random point cloud (should always work so no need to test seed)
  cloud.data <- cbind(runif(50), runif(50))
  
  # create equivalent distance matrix
  matrix.data <- calc.distmat(cloud.data)
  
  # get persistent homology for both
  phom.cloud <- calculate_homology(cloud.data)
  phom.matrix<- calculate_homology(matrix.data, format = "distmat")
  
  # make sure both have same persistent homology
  expect_equal(phom.cloud, phom.matrix)
})

# test output format works as expected
test_that("Matrix and data frame outputs are equal", {
  # generate required phoms (both types)
  data("unif2d")
  phom.mat <- calculate_homology(unif2d)
  phom.df <- calculate_homology(unif2d, return_df = TRUE)
  
  # compare
  expect_identical(as.matrix(phom.df), phom.mat)
})