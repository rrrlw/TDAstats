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
#' @param dataset numeric matrix containing point cloud or distance matrix
#' @param dim maximum dimension of features to calculate
#' @param threshold maximum diameter for computation of Vietoris-Rips complexes
#' @param type choose between "vr" for Vietoris-Rips complex or "cub" for
#'   cubical complex
#' @param p number of the prime field Z/pZ to compute the homology over
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
# EVENTUALLY NEED TO ADD ALPHA COMPLEX?
calculate_homology <- function(dataset, dim = 1, threshold = -1, type = "vr",
                               p = 2L, format = "cloud",
                               standardize = FALSE, return_df = FALSE) {
  
  # pick return format based on ripserr function params
  return_format <- ifelse(return_df, "df", "mat")
 
  # Deprecated: use ripserr::vietoris_rips or ripserr::cubical
  switch(type,
         "vr" = {
           .Deprecated(new = "ripserr::vietoris_rips")
           
           return(ripserr::vietoris_rips(dataset = dataset,
                                         dim = dim,
                                         threshold = threshold,
                                         p = p,
                                         format = format,
                                         standardize = standardize,
                                         return_format = return_format))
         },
         "cub" = {
           .Deprecated(new = "ripserr::cubical")
           
           return(ripserr::cubical(dataset = dataset,
                                   threshold = threshold,
                                   standardize = standardize,
                                   return_format = return_format))
         },
         stop(paste("Invalid value for parameter type (should be either `vr`",
         "or `cub`):", type)))
}
