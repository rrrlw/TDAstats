#####WASSERSTEIN CALCULATION#####
# convenience function using `wass_mat_calc` as main (`wass_workhorse` indirectly)
# as main workhorse
# ellipsis includes dim, format, threshold for calculate_homology
wass_cloud_calc <- function(pts1, pts2, pow.val = 1, ...) {
  # make sure pts1 and pts2 (matrices) have same # of cols
  if (ncol(pts1) != ncol(pts2)) {
    stop("Something wrong in code here; invalid arguments (unequal number of columns in point clouds.")
  }

  # calculate persistent homology
  phom1 <- calculate_homology(pts1, ...)
  phom2 <- calculate_homology(pts2, ...)
  
  # check if dim is present as parameter within ...
  dim <- 1 # default
  list.version <- list(...)
  if (length(list.version) > 0) {
    if (!is.null(list.version$dim)) {
      dim <- as.integer(list.version$dim)
    }
  }

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

# get top n features from each dimension
phom.topn <- function(phom, n) {
  # all the dimensions
  dims <- unique(phom[, 1])
  phom.adjust <- matrix(NA, nrow = 0, ncol = 3)
  
  # go through all dimensions
  for (curr.dim in dims) {
    # get top n rows from this dimension
    curr.phom <- phom[phom[, 1] == curr.dim, ]
    names(curr.phom) <- NULL
    curr.phom <- matrix(curr.phom, ncol = 3)
    
    # order
    curr.phom <- cbind(curr.phom,
                       curr.phom[, 3] - curr.phom[, 2])
    curr.phom <- curr.phom[order(curr.phom[, 4], decreasing = TRUE), ]
    curr.phom <- matrix(curr.phom, ncol = 4)
    
    # combine phom.adjust - guaranteed at least one row
    phom.adjust <- rbind(phom.adjust,
                         curr.phom[1:min(nrow(curr.phom), n), 1:3])
  }
  
  return(phom.adjust)
}

#####BOOTSTRAPPING#####
#' Identify Significant Features in Persistent Homology
#' 
#' An empirical method (bootstrap) to differentiate between features that
#' constitute signal versus noise based on the magnitude of their
#' persistence relative to one another. Note: you must have at
#' least 5 features of a given dimension to use this function.
#' 
#' @param features 3xn data frame of features; the first column must be
#'   dimension, the second birth, and the third death
#' @param dim dimension of features of interest
#' @param reps  number of replicates
#' @param cutoff  percentile cutoff past which features are considered
#'   significant
#' @examples
#' # get dataset (noisy circle) and calculate persistent homology
#' angles <- runif(100, 0, 2 * pi)
#' x <- cos(angles) + rnorm(100, mean = 0, sd = 0.1)
#' y <- sin(angles) + rnorm(100, mean = 0, sd = 0.1)
#' annulus <- cbind(x, y)
#' phom <- calculate_homology(annulus)
#' 
#' # find threshold of significance
#' # expecting 1 significant feature of dimension 1 (Betti-1 = 1 for annulus)
#' thresh <- id_significant(features = as.data.frame(phom),
#'                          dim = 1,
#'                          reps = 500,
#'                          cutoff = 0.975)
#'
#' # generate flat persistence diagram
#' # every feature higher than `thresh` is significant
#' plot_persist(phom, flat = TRUE)
#' @export
id_significant <- function(features, dim = 1,
                           reps = 100, cutoff = 0.975) {
  # subset only features of dimension of interest
  colnames(features) <- c("dimension", "birth", "death")
  features <- features[features[, 1] == dim, ]
  
  # make sure conditions are met
  if (nrow(features) < 3) {
    stop(paste("There are too few (< 3) features of the",
               "desired dimension in the feature matrix."))
  }
  
  # calculate persistence
  features$persist <- features$death - features$birth
  
  # do bootstrap
  ans <- numeric(reps)
  for (i in 1:reps) {
    curr_sample <- sample(x = features$persist,
                          size = nrow(features),
                          replace = TRUE)
    ans[i] <- mean(curr_sample)
  }
  
  # return threshold above which features count as "significant"
  stats::quantile(ans, cutoff, names = FALSE)
}

#####DISTANCE BETWEEN PERSISTENT HOMOLOGY#####
# calculate distance between two persistent homology matrices (filled with features)
# inspired from Wasserstein/EMD but not the same (no intrinsic ordering of features, etc.)
#' Calculate Distance between Homology Matrices
#' 
#' Calculates the distance between two matrices containing persistent homology
#' features, usually as returned by the `calculate_homology` function.
#' 
#' Note that the absolute value of this measure of distance is not meaningful
#' without a null distribution or at least another value for relative
#' comparison (e.g. finding most similar pair within a triplet).
#' 
#' @param phom1 3-by-n numeric matrix containing persistent homology for first dataset
#' @param phom2 3-by-n numeric matrix containing persistent homology for second dataset
#' @param limit.num limit comparison to only top `limit.num` features in each dimension
#' @return distance vector (1 element per dimension) between `phom1` and `phom2`
#' @export
phom.dist <- function(phom1, phom2, limit.num = 0) {
  # make sure both matrices have at least some features
  if (nrow(phom1) < 1 | nrow(phom2) < 1) {
    stop("Each homology matrix must have at least one feature.")
  }
  
  # make sure homology matrix format is correct
  if (ncol(phom1) != 3 | ncol(phom2) != 3) {
    stop("Each homology matrix must have exactly three columns.")
  }
  
  # none of the values in the homology matrix should be negative
  if (sum(phom1 < 0) > 0 |
      sum(phom2 < 0) > 0) {
    stop("A homology matrix cannot contain any negative values.")
  }
  
  # make sure limit.num is an integer
  if (!(class(limit.num) %in% c("integer", "numeric"))) {
    stop(limit.num, " must be of class integer or numeric.")
  }
  limit.num <- as.integer(limit.num)
  
  # select only top n features from each dimension if `limit.num` parameter > 0
  if (limit.num > 0) {
    phom1 <- phom.topn(phom1, n = limit.num)
    phom2 <- phom.topn(phom2, n = limit.num)
  }
  
  # calculate maximum feature dimension for each matrix
  max.dim1 <- max(phom1[, 1])
  max.dim2 <- max(phom2[, 1])

  # call function that does actual work
  wass_mat_calc(phom1, phom2, pow.val = 1, dim = max(max.dim1, max.dim2))
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
#' @param update if greater than zero, will print a message every `update` iterations
#' @param ... arguments for `calculate_homology` used for each permutation; this includes the `format`, `dim`, and `threshold` parameters
#' @return list containing results of permutation test
#' @export
permutation_test <- function(data1, data2, iterations,
                             exponent = 1, update = 0, ...) {
  # doesn't work for distance matrices, only point clouds; set dim parameter
  dim <- 1
  list.version <- list(...)
  if (length(list.version) > 0) {
    if (!is.null(list.version$format)) {
      if (list.version$format != "cloud") {
        stop("Permutation tests only work for point clouds.")
      }
    }
    # set dim parameter from ... that user passed
    if (!is.null(list.version$dim)) {
      dim <- as.integer(list.version$dim)
    }
  }
  
  if (update > 0) {
    cat("Starting function\n")
  }
  
  # make sure both are matrices with same number of columns,
  # sufficient number of rows, and no missing values
  if (!("matrix" %in% class(data1)) |
      !("matrix" %in% class(data2))) {
    stop("Both point clouds must be passed as matrices.")
  }
  if (update > 0) {
    cat("PASSED error check #1\n")
  }
  if (nrow(data1) < 2 | nrow(data2) < 2) {
    stop("Both point clouds must have at least 2 points (rows) each.")
  }
  if (update > 0) {
    cat("PASSED error check #2\n")
  }
  if (ncol(data1) < 2 | ncol(data2) < 2) {
    stop("Both point clouds must be in least 2 dimensions (columns) each.")
  }
  if (update > 0) {
    cat("PASSED error check #3\n")
  }
  if (ncol(data1) != ncol(data2)) {
    stop("Both point clouds must have the same number of dimensions.")
  }
  if (update > 0) {
    cat("PASSED error check #4\n")
  }
  if (sum(is.na(data1)) > 0 | sum(is.na(data2)) > 0) {
    stop("There should be no NAs in the point clouds passed to this function.")
  }
  if (update > 0) {
    cat("PASSED error check #5\n")
  }
  class.data1 <- class(data1[1,1])
  class.data2 <- class(data2[1,1])
  allowed.classes <- c("numeric", "integer")
  if (!(class.data1 %in% allowed.classes) |
      !(class.data2 %in% allowed.classes)) {
    stop("Point clouds must be formatted as matrices filled with integers or numerics.")
  }
  if (update > 0) {
    cat("PASSED error check #6\n")
  }
  if (iterations <= 1) {
    stop("Permutation test must have at least 2 iterations (preferably more).")
  }
  if (update > 0) {
    cat("PASSED error check #7\n")
  }
  if (update > 0) {
    cat("Beginning calculations\n")
  }

  # calculate Wasserstein values for actual point clouds (prior to permuting)
  orig.wass <- wass_cloud_calc(data1, data2, exponent, ...)
  
  if (update > 0) {
    cat("Initial calculation complete; starting iterations\n")
  }

  # calculate Wasserstein values for each permutation
  combo.pts <- rbind(data1, data2)
  wass.values <- matrix(-1, ncol = ncol(data1), nrow = iterations)
  worked <- vapply(X = 1:iterations,
                   FUN.VALUE = logical(1),
                   FUN = function(curr.iter) {
                     # print update message if necessary
                     if (update > 0 & curr.iter %% update == 0) {
                       cat("At iteration #", curr.iter, "\n", sep = "")
                     }
                     
                     # permute the point clouds
                     curr.pts <- combo.pts[sample.int(n = nrow(combo.pts),
                                                      size = nrow(combo.pts),
                                                      replace = FALSE), ]

                     # calculate Wasserstein for permuted clouds
                     curr.wass.calc <- wass_cloud_calc(pts1 = curr.pts[1:nrow(data1), ],
                                                       pts2 = curr.pts[(nrow(data1) + 1):nrow(curr.pts), ],
                                                       pow.val = exponent,
                                                       ...)
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
  if (update > 0) {
    cat("Completed calculations\n")
  }
  return(answer)
}
