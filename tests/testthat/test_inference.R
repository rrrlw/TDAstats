context("Test statistical inference on persistent homology")
library(TDAstats)

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
  
  # make sure it's commutative
  perm.test <- permutation_test(unif2d, circle2d, iterations = N.ITER, update = 0)
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

test_that("permutation_test detects invalid parameters correctly", {
  # reproducibility not needed cuz exactly values aren't being tested
  # incorrect format (lower distance matrix instead of point cloud)
  cloud.1 <- cbind(runif(50), runif(50))
  cloud.2 <- cbind(rnorm(50), rnorm(50))
  expect_error(permutation_test(cloud.1, cloud.2, iterations = 10, format = "distmat"),
               "Permutation tests only work for point clouds\\.")
  
  # error because passing data frame instead of matrix
  cloud.df <- as.data.frame(cloud.1)
  expect_error(permutation_test(cloud.df, cloud.2, iterations = 10),
               "Both point clouds must be passed as matrices\\.")
  
  # error because not enough points
  cloud.temp <- matrix(cloud.1[1, ], nrow = 1)
  expect_error(permutation_test(cloud.temp, cloud.2, iterations = 10),
               "Both point clouds must have at least 2 points \\(rows\\) each\\.")
  
  # error because not enough dimensions
  cloud.temp <- matrix(cloud.2[, 1], ncol = 1)
  expect_error(permutation_test(cloud.1, cloud.temp, iterations = 10),
               "Both point clouds must be in least 2 dimensions \\(columns\\) each\\.")
  
  # error because of dimension mismatch
  cloud.1 <- cbind(cloud.1, runif(50))
  expect_error(permutation_test(cloud.1, cloud.2, iterations = 10),
               "Both point clouds must have the same number of dimensions\\.")
  cloud.1 <- cloud.1[, 1:2]
  
  # error because of missing values
  temp <- as.numeric(cloud.1[1, 1]) # avoid pass-by-reference just in case (make actual copy)
  cloud.1[1, 1] <- NA
  expect_error(permutation_test(cloud.1, cloud.2, iterations = 10),
               "There should be no NAs in the point clouds passed to this function\\.")
  cloud.1[1, 1] <- temp
  
  # error because matrix does not contain only numeric/integer
  temp <- as.numeric(cloud.2[5, 2])
  cloud.2[5, 2] <- "hello"
  expect_error(permutation_test(cloud.1, cloud.2, iterations = 10),
               "Point clouds must be formatted as matrices filled with integers or numerics\\.")
  cloud.2 <- cbind(runif(50), runif(50))
  
  # error because insufficient iterations
  expect_error(permutation_test(cloud.1, cloud.2, iterations = 1),
               "Permutation test must have at least 2 iterations \\(preferably more\\)\\.")
})

test_that("Bootstrap for identification of significant features is working", {
  # do bootstrap on an annulus
  angles <- runif(100, 0, 2 * pi)
  x <- cos(angles) + rnorm(100, 0, 0.1)
  y <- sin(angles) + rnorm(100, 0, 0.1)
  annulus <- cbind(x, y)
  
  phom <- calculate_homology(annulus, return_df = TRUE)
  
  thresh <- id_significant(phom,
                           dim = 1,
                           reps = 500,
                           cutoff = 0.975)
  
  # check if exactly one feature is above the threshold
  phom <- phom[phom$dim == 1, ]
  phom$persist <- phom$death - phom$birth
  test_sol <- sum(phom$persist >= thresh)
  
  expect_equal(test_sol, 1)
})

test_that("Distance between persistent homology is calculated correctly", {
  # equal persistent homologies should have no difference
  phom.1 <- matrix(c(0, 0, 1,
                     1, 0, 1), ncol = 3, byrow = TRUE)
  phom.2 <- matrix(c(0, 0, 1,
                     1, 0, 1), ncol = 3, byrow = TRUE)
  ans <- phom.dist(phom.1, phom.2)
  names(ans) <- NULL  # remove names b/c of testthat::expect_equal check
  expect_equal(ans[1], 0)
  expect_equal(ans[2], 0)
  
  # one unit difference should have phom.dist equal to 1 (basic non-zero test)
  # also one dimension shouldn't affect result of another
  phom.2[2, 3] <- 2
  ans <- phom.dist(phom.1, phom.2)
  names(ans) <- NULL  # remove names b/c of testthat::expect_equal check
  expect_equal(ans[1], 0)
  expect_equal(ans[2], 1)
  
  # now both dimensions have non-zero values
  phom.1[1, 3] <- 3
  ans <- phom.dist(phom.1, phom.2)
  names(ans) <- NULL  # remove names b/c of testthat::expect_equal check
  expect_equal(ans[1], 2)
  expect_equal(ans[2], 1)
  
  # should throw error (too few rows)
  phom.temp <- matrix(1, nrow = 0, ncol = 3)
  expect_error(phom.dist(phom.temp, phom.2),
               "Each homology matrix must have at least one feature\\.")
  
  # should throw error (wrong number of columns)
  phom.temp <- matrix(1, nrow = 2, ncol = 4)
  expect_error(phom.dist(phom.temp, phom.1),
               "Each homology matrix must have exactly three columns\\.")
  
  # should throw error (negative value in second parameter)
  phom.temp <- matrix(c(0, 0, 1,
                        1, -1, 0), ncol = 3, byrow = TRUE)
  expect_error(phom.dist(phom.1, phom.temp),
               "A homology matrix cannot contain any negative values\\.")
  
  # make sure phom.dist makes a difference (<= phom.dist without limit.num)
  data("unif2d")
  data("circle2d")
  phom.unif <- calculate_homology(unif2d, dim = 1)
  phom.circ <- calculate_homology(circle2d, dim = 1)
  
  dists <- phom.dist(phom.unif, phom.circ)
  dists2<- phom.dist(phom.unif, phom.circ, limit.num = 2)
  expect_lte(dists2[1], dists[1])
  expect_lte(dists2[2], dists[2])
})