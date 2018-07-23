## create a 2-d unit circle with 100 points uniformly distributed on its circumference.
## serves as an example of a 2-dimensional point cloud with a prominent 1-cycle.
# set seed for reproducibility
set.seed(1)

# create dataset
angles <- runif(100, 0, 2 * pi)
circle2d <- cbind(cos(angles), sin(angles))

# add to package using devtools
devtools::use_data(circle2d)
