# TDAstats: topological data analysis in R <img src="man/figures/HexTDA.png" align="right" height="175" width="151"/>

[![Travis-CI Build Status](https://travis-ci.org/rrrlw/TDAstats.svg?branch=master)](https://travis-ci.org/rrrlw/TDAstats)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rrrlw/TDAstats?branch=master&svg=true)](https://ci.appveyor.com/project/rrrlw/TDAstats)
[![Coverage Status](https://img.shields.io/codecov/c/github/rrrlw/TDAstats/master.svg)](https://codecov.io/github/rrrlw/TDAstats?branch=master)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN version](http://www.r-pkg.org/badges/version/TDAstats)](https://CRAN.R-project.org/package=TDAstats)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/TDAstats)](https://CRAN.R-project.org/package=TDAstats)

[![JOSS DOI](http://joss.theoj.org/papers/10.21105/joss.00860/status.svg)](https://doi.org/10.21105/joss.00860)
[![Zenodo DOI](https://zenodo.org/badge/130141540.svg)](https://zenodo.org/badge/latestdoi/130141540)

## Overview

TDAstats is an R pipeline for computing persistent homology in topological data analysis.

## Installation

To install TDAstats, run the following R code:
```r
# install from CRAN
install.packages("TDAstats")

# install development version from GitHub
devtools::install_github("rrrlw/TDAstats")

# install development version with vignettes/tutorials
devtools::install_github("rrrlw/TDAstats", build_vignettes = TRUE)
```

## Sample code

The following sample code creates two synthetic datasets, and calculates and visualizes their persistent homology to showcase the use of TDAstats.

```r
# load TDAstats
library("TDAstats")

# load sample datasets
data("unif2d")
data("circle2d")

# calculate persistent homology for both datasets
unif.phom <- calculate_homology(unif2d, dim = 1)
circ.phom <- calculate_homology(circle2d, dim = 1)

# visualize first dataset as persistence diagram
plot_persist(unif.phom)

# visualize second dataset as topological barcode
plot_barcode(circ.phom)
```

A more detailed tutorial can be found in the package vignettes or at [this Gist](https://gist.github.com/rrrlw/2fd22a834a883cb66454b1dabab9fdcb).

## Functionality

TDAstats has 3 primary goals:

1.  *Calculation of persistent homology*: the C++
[Ripser](https://github.com/Ripser/ripser)
project is a lightweight library for calculating persistent homology
that outpaces all of its competitors. Given the importance of computational
efficiency, TDAstats naturally uses Ripser behind the scenes for homology
calculations, using the Rcpp package to integrate the C++ code into an R
pipeline (Ripser for R).

2.  *Statistical inference of persistent homology*: persistent homology can be
used in hypothesis testing to compare the topological structure of two point
clouds. TDAstats uses a permutation test in conjunction with the Wasserstein
metric for nonparametric statistical inference.

3.  *Visualization of persistent homology*: persistent homology is visualized
using two types of plots - persistence diagrams and topological barcodes.
TDAstats provides implementations of both plot types using the
[ggplot2](https://github.com/tidyverse/ggplot2)
framework. Having ggplot2 underlying the plots confers many advantages to the
user, including generation of publication-quality plots and customization using
the ggplot object returned by TDAstats.

## Contribute

To contribute to TDAstats, you can create issues for any bugs/suggestions on the [issues page](https://github.com/rrrlw/TDAstats/issues). You can also fork the TDAstats repository and create pull requests to add features you think will be useful for users.

## Applications and use cases

TDAstats mentions in the wild:

* M\'{e}moli F, Singhal K. A Primer on Persistent Homology of Finite Metric Spaces. Bulletin of Mathematical Biology. 2019; 81(7): 2074. Links to [paper](https://link.springer.com/article/10.1007/s11538-019-00614-z) and [preprint](https://arxiv.org/abs/1905.13400)
* Srinivasan R, Chander A. Understanding Bias in Datasets using Topological Data Analysis. Fujitsu Laboratories of America. 2019. [Link](http://ceur-ws.org/Vol-2419/paper_9.pdf)
* Rickert J. A Mathematician's Perspective on Topological Data Analysis and R. 2018. [Link](https://rviews.rstudio.com/2018/11/14/a-mathematician-s-perspective-on-topological-data-analysis-and-r/)
* [Blog post on Data Management](https://www.kaisataipale.net/blog/2019/02/22/data-management/)
* [Analyzing finance data](https://github.com/kaitai/Example-with-TDAstats)
* [R package for visualizing persistent homology](https://github.com/rrrlw/ggtda)
