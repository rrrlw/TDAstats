---
title: 'TDAstats: R pipeline for computing persistent homology in topological data analysis'
tags:
  - R
  - topological data analysis
  - persistent homology
  - Vietoris-Rips complex
  - statistical resampling
authors:
  - name: "Raoul R. Wadhwa"
    orcid: 0000-0003-0503-9580
    affiliation: "1, 2"
  - name: "Drew F.K. Williamson"
    orcid: 0000-0003-1745-8846
    affiliation: "1, 3"
  - name: "Andrew Dhawan"
  	orcid: 0000-0002-5027-1277
  	affiliation: "1, 4"
  - name: "Jacob G. Scott"
  	orcid: 0000-0003-2971-7673
  	affiliation: "1, 2, 3"
affiliations:
 - name: Cleveland Clinic Lerner College of Medicine, Case Western Reserve University, Cleveland, OH 44195, USA
   index: 2
 - name: Case Western Reserve University School of Medicine, Cleveland, OH 44106, USA
   index: 3
 - name: Department of Translational Hematology and Oncology Research, Cleveland Clinic Foundation, Cleveland, OH 44195, USA
   index: 1
 - name: Neurological Institute, Cleveland Clinic Foundation, Cleveland, OH 44195, USA
   index: 4
date: 27 July 2018
bibliography: paper.bib
---

# Summary

High-dimensional datasets are becoming more common in a variety of scientific fields. Well-known examples include next-generation sequencing in biology, patient health status in medicine, and computer vision in deep learning. Dimension reduction, using methods like principal component analysis (PCA), is a common preprocessing step for such datasets. However, while dimension reduction can save computing and human resources, it comes with the cost of significant information loss. Topological data analysis (TDA) aims to analyze the "shape" of high-dimensional datasets, without dimension reduction, by extracting features that are robust to small perturbations in data. Persistent features of a dataset can be used to describe it, and to compare it to other datasets. Visualization of persistent features can be done using topological barcodes or persistence diagrams (see figure). Application of TDA methods has granted greater insight into high-dimensional data [@tda-immune]; one prominent example of this is its use to characterize a clinically relevant subgroup of breast cancer patients [@tda-cancer].

The ``TDAstats`` R package is a comprehensive pipeline for conducting TDA. Once data is loaded into R, ``TDAstats`` can calculate, visualize, and conduct nonparametric statistical inference on persistent homology. The Ripser C++ library [@Ripser], benchmarked at approximately 40 times faster than comparable software, is wrapped using Rcpp [@rcpp-paper] for efficient computation of persistent homology. ``TDAstats`` generates topological barcodes and persistence diagrams using the ubiquitous ggplot2 library [@ggplot2-book], allowing use of ggplot2 functions to manipulate plots. This reduces the number of manual steps required to prepare publication-quality figures, thus enabling reproducible research [@reproduce-res]. ``TDAstats`` also implements nonparametric hypothesis testing of persistent homology using a permutation test, first described by @tda-testing. To our knowledge, ``TDAstats`` is the first library to implement this feature.

The primary barrier to using TDA is not mathematical comprehension. Although the algebraic topology that underlies TDA requires graduate-level study, the concepts necessary for application of TDA are far more intuitive. Rather, the barrier to entry is the lack of accessible, user-friendly software. ``TDAstats`` has an easy-to-use API with only 4 functions, each with only one or two intuitive parameters. Additionally, the provided vignettes cover its functionality with a comprehensive introduction and case study. Thus, even minimal knowledge of R will be sufficient to conduct TDA. We intend to use ``TDAstats`` to improve digit recognition algorithms, and hope that with its efficient implementation and user-friendly API, a far larger set of students and researchers can now apply TDA to answer research questions.

[Topological barcode (left) and persistence diagram (right) of the sphere3d sample dataset included with ``TDAstats``. The 0-cycles are colored red, the 1-cycles are colored green, and the 2-cycles are colored blue. For details on interpreting these plots, see @roadmap-tda.](sphere3d.png)

# Acknowledgments

We thank Francesca Buffa, Ph.D. and Chad Topaz, Ph.D. for their advice.

# References
