## set of 100 points uniformly distributed within a unit cube.
## serves as an example of a 3-dimensional point cloud with no
## prominent 0-, 1-, or 2-cycles
# set seed for reproducibility
set.seed(1)

# create dataset
unif3d <- cbind(runif(100), runif(100), runif(100))

# add to package using devtools
devtools::use_data(unif3d)
