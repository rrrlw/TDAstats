#' Calculate Persistent Homology of a Point Cloud
#'
#' Calculates the persistent homology of a point cloud, as represented by
#' a Vietoris-Rips complex. This function is an R wrapper for Ulrich Bauer's
#' Ripser C++ library for calculating persistent homology. For more
#' information on the C++ library, see https://github.com/Ripser/ripser.
#'
#' The `mat` parameter should be a numeric matrix with each row corresponding
#' to a single point, and each column corresponding to a single dimension. Thus,
#' if `mat` has 50 rows and 5 columns, it represents a point cloud with 50 points
#' in 5 dimensions.
#'
#' @param mat numeric matrix containing point cloud
#' @return 3-column matrix, with each row representing a TDA feature
#' @importFrom stats complete.cases
#' @export
#' @examples
#'
#' # create a 2-d point cloud of a circle (100 points)
#' num.pts <- 100
#' rand.angle <- runif(num.pts, 0, 2*pi)
#' pt.cloud <- cbind(cos(rand.angle), sin(rand.angle))
#'
#' # calculate persistent homology (num.pts by 3 numeric matrix)
#' pers.hom <- calculate_homology(pt.cloud)
calculate_homology <- function(mat) {

  # make sure matrix has at least 2 columns and at least 2 rows
  if (nrow(mat) < 2 | ncol(mat) < 2) {
    stop("Point cloud must have at least 2 points and at least 2 dimensions.")
  }

  # make sure matrix contains numeric (or integer) values
  # assumption: matrix can only hold objects of one class so only need to check
  #   one element
  temp <- mat[1, 1]
  if (class(temp) != "numeric" && class(temp) != "integer") {
    stop("Point cloud must contain values of class `numeric` or `integer` only.")
  }

  # make sure there are no NAs in matrix
  if (sum(stats::complete.cases(mat)) < nrow(mat)) {
    stop("Point cloud has missing values.")
  }

  # actually do work
  ans_vec <- ripser_cpp(mat)

  # format properly and return
  ans_mat <- matrix(ans_vec,
                    byrow = TRUE,
                    ncol = 3)
  colnames(ans_mat) <- c("dimension", "birth", "death")
  return(ans_mat)
}
