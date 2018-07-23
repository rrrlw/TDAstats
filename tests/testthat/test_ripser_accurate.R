context("Homology detection of significant features")
library(TDAstats)

test_that("calculate_homology detects 0- and 1-cycles in circle", {
  # don't need reproducibility here; with enough points,
  # it should always detect at least one 1-cycle (at least one 0-cycle always)
  
  # make dataset (2-d circle)
  angles <- runif(100, 0, 2 * pi)
  circle.data <- cbind(cos(angles), sin(angles))
  
  # calculate homology
  circle.hom <- calculate_homology(circle.data, dim = 1)
  
  # make sure at least one 0- and 1-cycle is detected (and only those)
  vals <- unique(circle.hom[, "dimension"])
  expect_equal(2, length(vals))
  expect_equal(0, min(vals))
  expect_equal(1, max(vals))
})

test_that("calculate_homology detects 0-, 1-, and 2-cycles in sphere", {
  # not worth running if 2-d test above works, no point wasting CRAN/Travis resources
  skip_on_cran()
  skip_on_travis()
  
  # reproducibility not needed (sufficient points should lead to 2-cycle)
  
  # make dataset (3-d sphere w/ Marsaglia method)
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
  
  # calculate homology for sphere
  sphere.hom <- calculate_homology(sphere.data, dim = 2)
  
  # make sure at least one 0-, 1-, and 2-cycle is detected (and only those)
  vals <- unique(sphere.hom[, "dimension"])
  #print(tail(sphere.hom))
  #print("---")
  #print(vals)
  expect_equal(3, length(vals))
  expect_equal(0, min(vals))
  expect_equal(1, sort(vals)[2])
  expect_equal(2, max(vals))
})