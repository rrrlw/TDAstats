# make sure matrix is nx3 (first col integer, next two numeric or integer)
#   and contains valid features (birth before death)
validate_matrix <- function(feature.matrix) {
  # exactly 3 columns: dimension, birth, death
  if (ncol(feature.matrix) != 3) {
    stop(paste("Invalid feature matrix with", ncol(feature.matrix),
               "columns: must have exactly 3 columns."))
  }
  # at least 1 row (at least 1 feature)
  if (nrow(feature.matrix) < 1) {
    stop("Invalid feature matrix: must have at least one row.")
  }
  # no NAs
  count.na <- sum(is.na(feature.matrix))
  if (count.na > 0) {
    stop(paste("Invalid feature matrix with", count.na,
               "NAs: should have zero NAs."))
  }
  # first column should be of type integer
  suppressWarnings(temp <- as.integer(feature.matrix[, 1]))
  if (sum(is.na(temp)) > 0) {
    stop("Invalid feature matrix: first column must be of type integer to represent feature dimension.")
  }
  # second and third columns should be of type numeric
  suppressWarnings(temp <- as.numeric(feature.matrix[, 2:3]))
  if (sum(is.na(temp)) > 0) {
    stop("Invalid feature matrix: second and third columns must be of type numeric to represent feature birth and death.")
  }
  # make sure all births are before corresponding deaths
  check.before <- sum(feature.matrix[, 2] > feature.matrix[, 3])
  if (check.before > 0) {
    stop("Invalid feature matrix: birth must come before death for all features.")
  }

  # valid matrix
  return(TRUE)
}
