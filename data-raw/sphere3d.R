## create a 3-d unit sphere with 100 points uniformly distributed on its surface.
## serves as an example of a 3-dimensional point cloud with a prominent 2-cycle.
# set seed for reproducibility
set.seed(1)

# Marsaglia method for picking points uniformly from sphere surface
rsphere <- function(num.points, radius = 1) {
  # generate points
  ans <- matrix(NA, nrow = num.points, ncol = 3)
  tmp <- vapply(1:num.points,
                FUN.VALUE = logical(1),
                FUN = function(curr.row) {
                  # pick appropriate x1 and x2
                  x1 <- 1
                  x2 <- 1
                  while (x1 * x1 + x2 * x2 >= 1) {
                    x1 <- runif(1, -1, 1)
                    x2 <- runif(1, -1, 1)
                  }
                  
                  # calculate values
                  x.val <- 2 * x1 * sqrt(1 - x1 ^ 2 - x2 ^ 2)
                  y.val <- 2 * x2 * sqrt(1 - x1 ^ 2 - x2 ^ 2)
                  z.val <- 1 - 2 * (x1 ^ 2 + x2 ^ 2)
                  
                  # store into answer matrix
                  ans[curr.row, 1] <<- x.val * radius
                  ans[curr.row, 2] <<- y.val * radius
                  ans[curr.row, 3] <<- z.val * radius
                  
                  # exit
                  return(TRUE)
                })
  
  # return answer
  return(ans)
}

# create dataset
sphere3d <- rsphere(100)

# add to package using devtools
devtools::use_data(sphere3d)
