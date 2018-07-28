## JOSS reviewer checklist

* Requirements for paper.md: **all of the components below are included in paper.md**
	- list of authors of software
	- author affiliations
	- summary for non-specialist audience
	- statement of need/purpose of software
	- ongoing research projects using software
	- list of key references
* Software license: **plaintext LICENSE file present for GPL-3 (OSI-approved)**
* High-level documentation: **README.md file present with purpose, installation instructions (dependencies in DESCRIPTION file), and example usage (more detailed examples in vignettes folder)**
* API documentation: **all public functions have appropriate documentation; non-plotting functions have examples**
* Community guidelines
	- contribute to software: **included in README.md**
	- report issues: **included in DESCRIPTION and README.md (GitHub issues page)**
	- seek support: **email included in DESCRIPTION file**
* Functionality: **vignettes (compile Rmarkdown from vignettes folder, or use devtools::build_vignettes() after installing package. Upon request, can also email HTML of built vignettes) introduce functionality**
* Tests: **tests folder contains automated testing suite using testthat R package**
* Novelty of software: **TDA package in R does have overlapping functionality, but TDAstats calculates homology faster (see https://github.com/Ripser/ripser for details), allows modifications of plots through ggplot2 (more customizable, reproducible graph), and allows for hypothesis testing with permutation tests (which no TDA software known to us implements)**
* Development environment: **the R programming language and the corresponding IDE (Rstudio) are open-source and free**
