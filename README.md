# TDAstats: topological data analysis in R <img src="man/figures/HexTDA.png" align="right" height="175" width="151"/>

[![Lifecycle:
superseded](https://img.shields.io/badge/lifecycle-superseded-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#superseded)

[![Travis-CI Build Status](https://travis-ci.org/rrrlw/TDAstats.svg?branch=master)](https://travis-ci.org/rrrlw/TDAstats)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rrrlw/TDAstats?branch=master&svg=true)](https://ci.appveyor.com/project/rrrlw/TDAstats)
[![Coverage Status](https://img.shields.io/codecov/c/github/rrrlw/TDAstats/master.svg)](https://codecov.io/github/rrrlw/TDAstats?branch=master)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN version](http://www.r-pkg.org/badges/version/TDAstats)](https://CRAN.R-project.org/package=TDAstats)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/TDAstats)](https://CRAN.R-project.org/package=TDAstats)

[![JOSS DOI](http://joss.theoj.org/papers/10.21105/joss.00860/status.svg)](https://doi.org/10.21105/joss.00860)
[![Zenodo DOI](https://zenodo.org/badge/130141540.svg)](https://zenodo.org/badge/latestdoi/130141540)

## Supersedure

TDAstats has transitioned to legacy status: It may be maintained to work with new versions of R, but no new features or refactoring are planned. For the essential functionality described below, we recommend the more targeted packages [ripserr](https://github.com/tdaverse/ripserr) to compute persistent homology, [phutil](https://github.com/tdaverse/phutil) to compute Kantorovich/Wasserstein distances, [inphr](https://github.com/tdaverse/inphr) for statistical inference, and [ggtda](https://github.com/tdaverse/ggtda) for visualization, as well as the original [TDA](https://github.com/compTAG/r-tda) package for a comprehensive standalone toolkit.

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

## Citation

If you use TDAstats, please consider citing the following (based on use):
* **General use of TDAstats**: Wadhwa RR, Williamson DFK, Dhawan A, Scott JG. TDAstats: R pipeline for computing persistent homology in topological data analysis. *Journal of Open Source Software*. 2018; 3(28): 860. doi: [10.21105/joss.00860](https://doi.org/10.21105/joss.00860)
* **TDAstats to calculate persistent homology (Ripser)**: Bauer U. Ripser: Efficient computation of Vietoris-Rips persistence barcodes. 2019; *arXiv*: 1908.02518.
* **TDAstats to perform statistical test**: Robinson A, Turner K. Hypothesis testing for topological data analysis. *J Appl Comput Topol*. 2017; 1: 241.

## Real-world applications, use cases, and mentions

* Stenseke J. Persistent homology and the shape of evolutionary games. Journal of Theoretical Biology. 2021; 531: 110903. Link to [paper](https://www.sciencedirect.com/science/article/pii/S0022519321003222).
* Torres-Espin A, Haefeli J, Ehsanian R, et al. Topological network analysis of patient similarity for precision management of acute blood pressure in spinal cord injury. eLife. 2021; 10: e68015. Link to [paper](https://elifesciences.org/articles/68015).
* Somasundaram E, Litzler A, Wadhwa R, Owen S, Scott J. Persistent homology of tumor CT scans is associated with survival in lung cancer. Medical Physics. 2021; 48(11): 7043-7051. Link to [paper](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1002/mp.15255) and [preprint](https://www.medrxiv.org/content/10.1101/2020.12.06.20244863v1).
* Richardson M, Verma R, Singhania A, Tabone O, Das M, Rodrigue M, Leissner P, Woltmann G, Cooper A, O'Garra A, Haldar P. Blood transcriptional phenotypes of progressive latent M. tuberculosis infection inform novel signatures that improve prediction of tuberculosis risk. Cell Reports Medicine. 2021. Link to [paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3815986).
* Perez-Moraga R, Fores-Martos J, Suay-Garcia B, Duval J-L, Falco A, Climent J. A COVID-19 Drug Repurposing Strategy through Quantitative Homological Similarities Using a Topological Data Analysis-Based Framework. Pharmaceutics. 2021; 13(4): 488. Link to [paper](https://www.mdpi.com/1999-4923/13/4/488).
* Kandanaarachchi S, Hyndman RJ. Leave-one-out kernel density estimates for outlier detection. Monash University. 2021. Link to [paper](https://www.monash.edu/business/ebs/research/publications/ebs/wp02-2021.pdf).
* Somasundaram EV, Brown SE, Litzler A, Scott JG, Wadhwa RR. Benchmarking R packages for calculation of persistent homology. R Journal. 2021; 13(1): 184-193. Link to [paper](https://journal.r-project.org/archive/2021/RJ-2021-033/RJ-2021-033.pdf).
* Brochard A, Blaszczyszyn B, Mallat S, Zhang S. Particle gradient descent model for point process generation. 2020. arXiv:2010.14928. Link to [preprint](https://arxiv.org/abs/2010.14928).
* Nguyen DQN, Xing L, Lin L. Community detection, pattern recognition, and hypergraph-based learning: approches using metric geometry and persistent homology. 2020. arXiv:2010.00435. Link to [preprint](https://arxiv.org/abs/2010.00435).
* Pinto GVF. Motivic constructions on graphs and networks with stability results. Doctoral Thesis: Universidade Estadual Paulista Rio Claro & Ohio State University. 2020. Link to [thesis](http://hdl.handle.net/11449/192494).
* Gommel M. A Machine Learning Exploration of Topological Data Analysis Applied to Low and High Dimensional fMRI Data. Doctoral Thesis: University of Iowa. 2019. doi: [10.17077/etd.005247](https://doi.org/10.17077/etd.005247). Link to [thesis](https://iro.uiowa.edu/discovery/fulldisplay/alma9983779398602771/01IOWA_INST:ResearchRepository?tags=scholar).
* Mémoli F, Singhal K. A Primer on Persistent Homology of Finite Metric Spaces. Bulletin of Mathematical Biology. 2019; 81(7): 2074. Links to [paper](https://link.springer.com/article/10.1007/s11538-019-00614-z) and [preprint](https://arxiv.org/abs/1905.13400)
* Srinivasan R, Chander A. Understanding Bias in Datasets using Topological Data Analysis. Fujitsu Laboratories of America. 2019. [Link](http://ceur-ws.org/Vol-2419/paper_9.pdf)
* Kough D, Neuzil M, Simpson C, Glover R. Analyzing State of the Union Addresses using Topology. University of St. Thomas. 2019. [Link](https://www.stthomas.edu/media/collegeofartsandsciences/mathematics/pdf/camsummer2019/CAMReport2019TDANeuzilKoughSimpson.pdf)
* Rickert J. A Mathematician's Perspective on Topological Data Analysis and R. 2018. [Link](https://rviews.rstudio.com/2018/11/14/a-mathematician-s-perspective-on-topological-data-analysis-and-r/)
* [Blog post on Data Management](https://www.kaisataipale.net/blog/2019/02/22/data-management/)
* [Analyzing finance data](https://github.com/kaitai/Example-with-TDAstats)
* [R package for visualizing persistent homology](https://github.com/rrrlw/ggtda)
