## set of 100 points uniformly distributed within a unit square.
## serves as an example of a 2-dimensional point cloud with no
## prominent 0- or 1-cycles
# set seed for reproducibility
set.seed(1)

# create dataset
unif2d <- cbind(runif(100), runif(100))

# add to package using devtools
devtools::use_data(unif2d)
