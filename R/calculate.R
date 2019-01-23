#' Calculate Persistent Homology of a Point Cloud
#'
#' Calculates the persistent homology of a point cloud, as represented by
#' a Vietoris-Rips complex. This function is an R wrapper for Ulrich Bauer's
#' Ripser C++ library for calculating persistent homology. For more
#' information on the C++ library, see <https://github.com/Ripser/ripser>.
#'
#' The `mat` parameter should be a numeric matrix with each row corresponding
#' to a single point, and each column corresponding to a single dimension. Thus,
#' if `mat` has 50 rows and 5 columns, it represents a point cloud with 50 points
#' in 5 dimensions. The `dim` parameter should be a positive integer.
#' Alternatively, the `mat` parameter could be a distance matrix (upper
#' triangular half is ignored); note: `format` should be specified as "ldm".
#'
#' @param mat numeric matrix containing point cloud or distance matrix
#' @param dim maximum dimension of features to calculate
#' @param threshold maximum diameter for computation of Vietoris-Rips complexes
#' @param format  format of `mat`, either "cloud" for point cloud or "distmat" for distance matrix
#' @param standardize boolean determining whether point cloud size should be standardized
#' @param return_df defaults to `FALSE`, returning a matrix;
#'   if `TRUE`, returns a data frame
#' @return 3-column matrix or data frame, with each row representing a TDA feature
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
calculate_homology <- function(mat, dim = 1, threshold = -1, format = "cloud",
                               standardize = FALSE, return_df = FALSE) {
  # coerce mat into matrix to work with class object such as dist class object
  mat <- as.matrix(mat)
  
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
  
  # make sure dim is an integer greater than or equal to zero
  if (as.integer(dim) != dim) {
    stop("dim parameter needs to be an integer")
  }
  if (dim < 0) {
    stop("dim cannot be negative")
  }
  
  # make sure threshold is of type numeric
  if (!(class(threshold) %in% c("numeric", "integer"))) {
    stop("threshold parameter must be of type numeric")
  }
  threshold <- as.numeric(threshold)
  
  # make sure format is either "cloud" or "distmat"
  if (!(format %in% c("cloud", "distmat"))) {
    stop("format parameter should be either \"cloud\" or \"distmat\"")
  }
  int.format <- ifelse(format == "cloud", 0, 1)

  # standardize if necessary
  # method: independently for each coordinate axis: (x - min) / (max - min)
  if (standardize) {
    for (i in 1:ncol(mat)) {
      min.val <- min(mat[, i])
      max.val <- max(mat[, i])
      
      # skip if only one unique value in this column (all identical)
      if (min.val == max.val) next
      
      mat[, i] <- (mat[, i] - min.val) / (max.val - min.val)
    }
  }
  
  # actually do work
  ans_vec <- ripser_cpp(mat, dim, threshold, int.format)

  # format properly and return
  ans_mat <- matrix(ans_vec,
                    byrow = TRUE,
                    ncol = 3)
  colnames(ans_mat) <- c("dimension", "birth", "death")
  
  # convert matrix to data frame - better for visualization, manipulation
  if (return_df) {
    ans_df <- as.data.frame(ans_mat)
    ans_df$dimension <- as.integer(ans_df$dimension)
    return(ans_df)
  } else {
    return(ans_mat)
  }
}
