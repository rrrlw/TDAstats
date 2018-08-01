#####WASSERSTEIN CALCULATION#####
# convenience function using `wass_mat_calc` as main (`wass_workhorse` indirectly)
# as main workhorse
wass_cloud_calc <- function(pts1, pts2, pow.val = 1, dim = 1, format = "cloud") {
  # make sure pts1 and pts2 (matrices) have same # of cols
  if (ncol(pts1) != ncol(pts2)) {
    stop("Something wrong in code here; invalid arguments (unequal number of columns in point clouds.")
  }

  # calculate persistent homology
  phom1 <- calculate_homology(pts1, dim = dim, format = format)
  phom2 <- calculate_homology(pts2, dim = dim, format = format)

  # return Wasserstein metric
  return(wass_mat_calc(phom1, phom2, pow.val = pow.val, dim = dim))
}

# convenience function using `wass_workhorse` as main workhorse
# mat1 and mat2 are feature matrices for which we want to calculate the
# Wasserstein metric
wass_mat_calc <- function(mat1, mat2, pow.val = 1, dim = 1) {
  # only want to calculate Wasserstein for overlapping dimensions
  min.dim <- 0
  max.dim <- dim
  
  # make workhorse function do all the work
  wass.vals <- vapply(X = min.dim:max.dim,
                      FUN.VALUE = numeric(1),
                      FUN = function(curr.dim) {
                        # subset matrices for current dimension only
                        curr.mat1 <- matrix(mat1[mat1[, 1] == curr.dim, ], ncol = 3)
                        curr.mat2 <- matrix(mat2[mat2[, 1] == curr.dim, ], ncol = 3)

                        #print(dim(curr.mat1))
                        #if (nrow(curr.mat1) == 0 |
                        #    nrow(curr.mat2) == 0) {
                        #  stop("Calculating Wasserstein for invalid dimension")
                        #}
                        # if there's no features for a given dimension in current permutation
                        # add a 
                        if (nrow(curr.mat1) == 0) {
                          curr.mat1 <- matrix(c(curr.dim, 0, 0), ncol = 3)
                        }
                        if (nrow(curr.mat2) == 0) {
                          curr.mat2 <- matrix(c(curr.dim, 0, 0), ncol = 3)
                        }

                        # calculate feature lengths from feature matrices
                        curr.feat1 <- curr.mat1[, 3] - curr.mat1[, 2]
                        curr.feat2 <- curr.mat2[, 3] - curr.mat2[, 2]

                        # return Wasserstein metric for current dimension
                        this.ans <- wass_workhorse(curr.feat1, curr.feat2, pow.val)
                        if (is.na(this.ans)) stop("Error: found an NA somewhere1")
                        #print(paste("This ans:", this.ans))
                        return(this.ans)
                      })

  # return vector of Wasserstein metrics (one per dimension)
  names(wass.vals) <- as.character(min.dim:max.dim)
  return(wass.vals)
}

# vec1 and vec2 are numeric vectors containing the duration (death - birth) of
#   all the features (of the same dimension) of their respective calculated
#   persistent homologies
wass_workhorse <- function(vec1, vec2, pow.val = 1) {
  # make sure vec1 and vec2 are of the same length
  if (length(vec1) < length(vec2)) {
    vec1 <- c(vec1, rep(0, length(vec2) - length(vec1)))
  } else if (length(vec2) < length(vec1)) {
    vec2 <- c(vec2, rep(0, length(vec1) - length(vec2)))
  }

  # make sure vec1 and vec2 are now of equal length
  # user should never see this
  if (length(vec1) != length(vec2)) stop("Something wrong with length-equaling code above (inference.R; wass_calc).")
  if (length(vec1) == 0) stop("Error with feature calculating code.")
  
  # sort both
  vec1 <- sort(vec1, decreasing = TRUE)
  vec2 <- sort(vec2, decreasing = TRUE)

  # make vector of absolute value of differences
  diff.vec <- abs(vec1 - vec2)

  # return answer (incorporating power value parameter)
  ans <- sum(diff.vec ^ pow.val)
  if (is.na(ans)) {
    stop("Error: found an NA somewhere2")
  }
  return(ans)
}

#####PERMUTATION TEST#####
#' Statistical Inference for Topological Data Analysis
#'
#' Conducts a permutation test for nonparametric statistical inference
#' of persistent homology in topological data analysis.
#'
#' The persistent homology of two point clouds are compared with the
#' Wasserstein metric (where Wasserstein-1 is also known as the Earth
#' Mover's Distance). However, the magnitude of the metric for a single pair
#' of point clouds is meaningless without a reference distribution. This
#' function uses a permutation test (permuting the points between the two
#' clouds) as a nonparametric hypothesis test for statistical inference.
#'
#' For more details on permutation tests for statistical inference in
#' topological data analysis, see Robinson A, Turner K. Hypothesis
#' testing for topological data analysis. J Appl Comput Topology. 2017;
#' 1(2): 241-261.<doi:10.1007/s41468-017-0008-7>
#'
#' @param data1 first dataset
#' @param data2 second dataset
#' @param iterations number of iterations for distribution in permutation test
#' @param exponent parameter `p` that returns Wasserstein-p metric
#' @param dim maximum dimension of cycles for which to compare homology
#' @param format  format of data, either "cloud" for point cloud or "distmat" for distance matrix
#' @return list containing results of permutation test
#' @export
permutation_test <- function(data1, data2, iterations,
                             exponent = 1, dim = 1,
                             format = "cloud") {
  # make sure both are matrices with same number of columns,
  # sufficient number of rows, and no missing values
  if (class(data1) != "matrix" |
      class(data2) != "matrix") {
    stop("Both point clouds must be passed as matrices.")
  }
  if (nrow(data1) < 2 | nrow(data2) < 2) {
    stop("Both point clouds must have at least 2 points (rows) each.")
  }
  if (ncol(data1) < 2 | ncol(data2) < 2) {
    stop("Both point clouds must be in least 2 dimensions (columns) each.")
  }
  if (ncol(data1) != ncol(data2)) {
    stop("Both point clouds must have the same number of dimensions.")
  }
  if (sum(is.na(data1)) > 0 | sum(is.na(data2)) > 0) {
    stop("There should be no NAs in the point clouds passed to this function.")
  }
  class.data1 <- class(data1[1,1])
  class.data2 <- class(data2[1,1])
  allowed.classes <- c("numeric", "integer")
  if (!(class.data1 %in% allowed.classes) |
      !(class.data2 %in% allowed.classes)) {
    stop("Point clouds must be formatted as matrices filled with integers or numerics.")
  }
  if (iterations <= 1) {
    stop("Permutation test must have at least 2 iterations (preferably more).")
  }

  # calculate Wasserstein values for actual point clouds (prior to permuting)
  orig.wass <- wass_cloud_calc(data1, data2, exponent, dim = dim, format)

  # calculate Wasserstein values for each permutation
  combo.pts <- rbind(data1, data2)
  wass.values <- matrix(-1, ncol = ncol(data1), nrow = iterations)
  worked <- vapply(X = 1:iterations,
                   FUN.VALUE = logical(1),
                   FUN = function(curr.iter) {
                     # permute the point clouds
                     curr.pts <- combo.pts[sample.int(n = nrow(combo.pts),
                                                      size = nrow(combo.pts),
                                                      replace = FALSE), ]

                     # calculate Wasserstein for permuted clouds
                     curr.wass.calc <- wass_cloud_calc(pts1 = curr.pts[1:nrow(data1), ],
                                                       pts2 = curr.pts[(nrow(data1) + 1):nrow(curr.pts), ],
                                                       pow.val = exponent,
                                                       dim = dim,
                                                       format = format)
                     #print(curr.wass.calc)

                     # store into matrix
                     wass.values[curr.iter, ] <<- curr.wass.calc[1:ncol(wass.values)]

                     # all went okay
                     return(TRUE)
                   })

  # return list of lists; each smaller list corresponds to a single dimension
  # each smaller list has: (1) p-value (not corrected multiple testing)
  #                        (2) all trajectories (to plot as histogram)
  #                        (3) Wasserstein original metric for this dimension
  #                        (4) Title (character string with dimension to avoid mismatch w/ indexes)
  answer <- lapply(X = 0:dim,
                   FUN = function(curr.dim) {
                     curr.ans <- list()
                     curr.ans$dimension <- curr.dim
                     curr.ans$permvals <- wass.values[, curr.dim + 1]
                     curr.ans$wasserstein <- orig.wass[curr.dim + 1]
                     names(curr.ans$wasserstein) <- NULL
                     curr.ans$pvalue <- mean(curr.ans$permvals > curr.ans$wasserstein)

                     return(curr.ans)
                   })
  return(answer)
}
