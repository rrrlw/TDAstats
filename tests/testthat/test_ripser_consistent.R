context("Homology calculation consistency")
library(TDAstats)

test_that("calculate_homology is consistent with standardization", {
  # for reproducibility
  set.seed(0)
  
  # make dataset (random uniform in unit square)
  square.data <- cbind(runif(50), runif(50))
  
  # test homology calculation for 0-cycles and 1-cycles
  expect_equal_to_reference(calculate_homology(square.data, dim = 1,
                                               standardize = TRUE),
                            file = "2dconsist-std")
})

test_that("calculate_homology is consistent in 2-d", {
  # for reproducibility
  set.seed(1)
  
  # make dataset (random uniform in unit square)
  square.data <- cbind(runif(50), runif(50))
  
  # test homology calculation for 0-cycles and 1-cycles
  expect_equal_to_reference(calculate_homology(square.data, dim = 1),
                            file = "2dconsist-1")
  
  # make dataset (circle in 2-d)
  angles <- runif(50, 0, 2 * pi)
  circle.data <- cbind(cos(angles), sin(angles))
  
  # test homology calculation for 0-cycles and 1-cycles
  expect_equal_to_reference(calculate_homology(square.data, dim = 1),
                            file = "2dconsist-2")
})

test_that("calculate_homology is consistent in 3-d", {
  # dim = 2 calculations take far longer, not worth checking on CRAN or Travis
  # if above test already works
  skip_on_cran()
  skip_on_travis()
  
  # for reproducibility
  set.seed(1)
  
  # make dataset (random uniform in unit cube)
  cube.data <- cbind(runif(50), runif(50), runif(50))
  
  # test homology calculation for 0-cycles, 1-cycles, and 2-cycles
  expect_equal_to_reference(calculate_homology(cube.data, dim = 2),
                            file = "3dconsist-1")
  
  # make dataset (sphere in 2-d with Marsaglia method)
  sphere.data <- matrix(NA, nrow = 100, ncol = 3)
  vapply(1:nrow(sphere.data),
         FUN.VALUE = logical(1),
         FUN = function(curr.row) {
           # pick x1 and x2 from unif dist until condition met
           x1 <- 1
           x2 <- 1
           temp.triangle <- x1 * x1 + x2 * x2
           while (temp.triangle >= 1) {
             x1 <- runif(1, -1, 1)
             x2 <- runif(1, -1, 1)
             
             temp.triangle <- x1 * x1 + x2 * x2
           }
           
           # calculate coordinates
           common.calc <- sqrt(1 - temp.triangle)
           x.val <- 2 * x1 * common.calc
           y.val <- 2 * x2 * common.calc
           z.val <- 1 - 2 * temp.triangle
           
           # store into matrix and exit
           sphere.data[curr.row, 1] <<- x.val
           sphere.data[curr.row, 2] <<- y.val
           sphere.data[curr.row, 3] <<- z.val
           return(TRUE)
         })
  
  # test homology calculation for 0-cycles, 1-cycles, and 2-cycles
  expect_equal_to_reference(calculate_homology(sphere.data, dim = 2),
                            file = "3dconsist-2")
})