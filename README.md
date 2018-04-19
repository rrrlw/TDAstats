TDAstats
===========================================
[![Travis-CI Build Status](https://travis-ci.org/rrrlw/TDAstats.svg?branch=master)](https://travis-ci.org/rrrlw/TDAstats)

Overview
--------

TDAstats is an R pipeline for topological data analysis, specifically, the
use of persistent homology in Vietoris-Rips simiplicial complexes to study the
shape of data.

Installation
------------

To install TDAstats, run the following R code:
```r
# install TDAstats from GitHub repository
devtools::install_github("rrrlw/TDAstats")
```
Note that you need the devtools package installed first. If you do not have
devtools installed, you can do so with the following R code:
```r
# install devtools from CRAN
install.packages("devtools")

# install TDAstats from GitHub repository
devtools::install_github("rrrlw/TDAstats")
```

We hope to submit TDAstats to CRAN soon. If accepted, the code above will be
revised to reflect installation from CRAN, rather than a GitHub repository.

Usage
-----

TDAstats has 3 primary goals:

1.  *Calculation of persistent homology*: the C++
[Ripser](https://github.com/Ripser/ripser)
project is a lightweight library for calculating persistent homology
that outpaces all of its competitors. Given the importance of computational
efficiency, TDAstats naturally uses Ripser behind the scenes for homology
calculations, using the Rcpp package to integrate the C++ code into an R
pipeline.

2.  *Statistical inference of persistent homology*: persistent homology can be
used in hypothesis testing to compare the topological structure of two point
clouds. TDAstats uses a permutation test in conjunction with the Wasserstein
metric for nonparametric statistical inference.

3.  *Visualization of persistent homology*: persistent homology is visualized
using two types of plots - persistence diagrams and topological barcodes.
TDAstats provides implementations of both plot types using the ggplot2
framework. Having ggplot2 underlying the plots confers many advantages to the
user, including generation of publication-quality plots and customization using
the ggplot object returned by TDAstats.
