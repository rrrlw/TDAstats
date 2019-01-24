#####PERSISTENCE DIAGRAM#####
# plot features as persistence diagram
# returns the plotted ggplot instance
#' Plot Persistent Homology as Persistence Diagram
#'
#' Plots a feature matrix as a persistence diagram. See `plot_barcode` for an
#' alternate visualization method of persistent homology.
#'
#' The `feature.matrix` parameter should be a numeric matrix with each
#' row corresponding to a single feature. It should have 3 columns
#' corresponding to feature dimension (col 1), feature birth (col 2), and
#' feature death (col 3). The first column should be filled with integers,
#' and the next two columns should be filled with numeric values.
#' The output from the `calculate_homology` function in this package will be a
#' valid value for the `feature.matrix` parameter.
#'
#' This function uses the ggplot2 framework to generate persistence diagrams.
#' For details, see: Wickham H (2009, ISBN:9780387981413). ggplot2: Elegant Graphics for Data
#' Analysis. Springer-Verlag: New York, NY.
#'
#' @param feature.matrix nx3 matrix representing persistent homology features
#' @param flat default FALSE; if TRUE, plots flat persistent homology instead
#' @param cutoff threshold for significant features; line added as marker on plot
#' @return ggplot instance representing persistence diagram
#' @import ggplot2
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
#'
#' # plot calculated homology features as persistence diagram
#' plot_persist(pers.hom)
plot_persist <- function(feature.matrix, flat = FALSE,
                         cutoff = 0) {
  # make sure feature matrix is formatted properly
  validate_matrix(feature.matrix)

  # get graphing parameters
  # N.B.: x- and y-axes have same min and max (always square graph)
  axes.min <- 0
  axes.max <- max(feature.matrix[, 3]) # don't have to check births

  # fix feature matrix formatting (into data frame for ggplot)
  names(feature.matrix) <- c("dimension", "birth", "death")
  feature.df <- as.data.frame(feature.matrix)
  feature.df$dimension <- as.factor(feature.df$dimension) # for colors to work correctly (discrete legend, not continuous color scale)

  # plot w/ ggplot
  if (!flat) {
    df.geompath <- data.frame(x = c(0, axes.max),
                              y = c(0, axes.max))
    g <- ggplot2::ggplot(data = feature.df) +
      ggplot2::xlim(axes.min, axes.max) + ggplot2::ylim(axes.min, axes.max) +                           # axis limits
      ggplot2::geom_abline(slope = 1, intercept = 0) +  # reference line
      #ggplot2::geom_path(data = df.geompath, ggplot2::aes_string(x = "x", y = "y")) +          # reference segment
      ggplot2::xlab("Feature appearance") + ggplot2::ylab("Feature disappearance") +                                 # axis titles
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),                               # add axis lines
                     panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),         # remove gridlines
                     panel.background = ggplot2::element_blank()) +                                           # remove default background color
     ggplot2::geom_point(ggplot2::aes_string(x = "birth", y = "death", shape = "dimension", colour = "dimension")) +   # add features as points
      ggplot2::coord_fixed(ratio = 1)
    
    # add line for bootstrap if desired
    if (cutoff > 0) {
      g +
        geom_abline(slope = 1, intercept = cutoff, linetype = 3, colour = "grey")
    } else {
      g
    }
  }
  else {
    feature.df$persistence <- feature.df$death - feature.df$birth
    g <- ggplot2::ggplot(data = feature.df) +
      ggplot2::geom_point(ggplot2::aes_string(x = "birth", y = "persistence", shape = "dimension", colour = "dimension")) +
      ggplot2::xlab("Feature appearance") +
      ggplot2::ylab("Feature persistence") +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                     panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    
    # add line for bootstrap if desired
    if (cutoff > 0) {
      g +
        geom_abline(slope = 0, intercept = cutoff, linetype = 3, colour = "grey")
    } else {
      g
    }
  }
}

#####TOPOLOGICAL BARCODE#####
# plot features as topological barcode
#' Plot Persistent Homology as Topological Barcode
#'
#' Plots a feature matrix as a topological barcode. See `plot_persist` for an
#' alternate visualization method of persistent homology.
#'
#' The `feature.matrix` parameter should be a numeric matrix with each
#' row corresponding to a single feature. It should have 3 columns
#' corresponding to feature dimension (col 1), feature birth (col 2), and
#' feature death (col 3). The first column should be filled with integers,
#' and the next two columns should be filled with numeric values.
#' The output from the `calculate_homology` function in this package will be a
#' valid value for the `feature.matrix` parameter.
#'
#' This function uses the ggplot2 framework to generate persistence diagrams.
#' For details, see: Wickham H (2009). ggplot2: Elegant Graphics for Data
#' Analysis. Springer-Verlag: New York, NY.
#'
#' @param feature.matrix nx3 matrix representing persistent homology features
#' @return ggplot instance representing topological barcode
#' @import ggplot2
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
#'
#' # plot calculated homology features as persistence diagram
#' plot_barcode(pers.hom)
plot_barcode <- function(feature.matrix) {
  # make sure feature matrix is formatted properly
  validate_matrix(feature.matrix)

  # get graphing parameters
  x.min <- 0
  x.max <- max(feature.matrix[, 3])
  y.min <- 1
  y.max <- nrow(feature.matrix)

  # fix formatting and add extra column as needed
  feature.df <- as.data.frame(feature.matrix)
  feature.df$dimension <- as.factor(feature.df$dimension) # for colors to work correctly (discrete legend, not continuous color scale)
  feature.df$vertical.pos <- y.min:y.max

  # plot w/ ggplot
  ggplot2::ggplot(data = feature.df) +
         ggplot2::xlim(x.min, x.max) + ggplot2::ylim(y.min, y.max) +                                                  # axis limits
         ggplot2::xlab("Vietoris-Rips Diameter") + ggplot2::ylab("") +                                                  # axis titles
         ggplot2::theme(axis.line.x = ggplot2::element_line(colour = "black"),                                        # add x-axis line
                        axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), # remove y-axis stuff
                        panel.grid = ggplot2::element_blank(),                                                        # remove gridlines
                        panel.background = ggplot2::element_blank()) +                                                # remove default background color
         ggplot2::geom_segment(ggplot2::aes_string(x = "birth", y = "vertical.pos",                                              # add actual bars for barcode
                          xend = "death", yend = "vertical.pos",
                          colour = "dimension"))
}
