---
title: 'Ripser.py: A Lean Persistent Homology Library for Python'
tags:
 - Python
 - Persistent Homology
 - Topological Data Analysis
authors:
 - name: Christopher Tralie
   orcid: 0000-0003-4206-1963
   affiliation: 1
 - name: Nathaniel Saul
   orcid: 0000-0002-8549-9810
   affiliation: 2
 - name: Rann Bar-On
   orcid: 0000-0002-4675-222X
   affiliation: 1
affiliations:
 - name: Department of Mathematics, Duke University
   index: 1
 - name: Department of Mathematics and Statistics, Washington State University
   index: 2
date: 10 September 2018
bibliography: paper.bib
---

# Summary
Topological data analysis (TDA) [@edelsbrunner2010computational], [@carlsson2009topology] is a field focused on understanding the shape and structure of data by computing topological descriptors that summarize features as connected components, loops, and voids. TDA has found wide applications across nonlinear time series analysis [@perea2015sliding], computer vision [@perea2014klein], computational neuroscience [@giusti2015clique], [@bendich2016persistent], computational biology [@iyer2010imaging], [@wu2017optimal], and materials science [@kramar2013persistence], to name a few of the many areas of impact in recent years.

Persistent homology [@edelsbrunner2010computational] is the main workhorse of TDA, and it computes a data structure known as the persistence diagram to summarize the space of stable topological features. The most commonly used scheme for generating persistence diagrams is the Vietoris Rips filtration (VR) since it is easily defined for any point cloud. In its naive implementation, VR is prohibitively slow, but recently a C++ library known as Ripser [@bauer2017ripser] has been devised to aggregate all known computational speedups of the VR filtration into one concise implementation. Because of the unprecedented speed of Ripser, it has created a large user base for both research and applications. However, the library as it stands is only a command line tool and as a result, multiple efforts have been made to wrap the C++ library for use in other languages, often via clunky system calls to the command line.

In this work, we develop an intuitive interface for VR filtrations with Ripser at its core via Cython. We have gone through extensive testing via continuous integration frameworks to ensure it works across all platforms and as a result, Ripser.py is currently as easy to setup as `pip install ripser`. We see this package as particularly useful for mathematicians with little programming experience who would like to use TDA as an entry point into data science, or conversely for researchers with little understanding of Algebraic Topology who would like to apply TDA to their problem domain. To aid this, we have created a large set of Jupyter notebooks to showcase some of the many applications that are possible with this library.

# Library Details
Ripser.py supplies two interfaces: one lightweight and functional interface, as well as an object-oriented interface designed to fit within the Scikit-Learn transformer paradigm [@scikit-learn]. We have merged together multiple branches of the original Ripser library ("sparse-distance-matrix," "representative-cocycles") to expose some lesser known but incredibly useful features hidden in Ripser.  Below we detail some of the special features made easy with our library.

## Sparse filtrations
Ripser.py can accomodate sparse distance matrices using the `scipy.sparse` library. In an accompanying notebooks, we demonstrate how Ripser.py can be used in conjunction with "sparse filtration" approximation algorithms [@cavannageometric] to add further computational speedups to Ripser.  Additionally, sparse matrices allow us to easily define a "sublevelset filtration," or the "watershed algorithm," for quantifying critical points in grid data such as time series and images, and Ripser.py includes a helper function to make this easy for image data. One of the notebooks demonstrates how this can be used to identify cells in an image, for instance.

## Coefficient Fields
Ripser.py naturally allows the specification of arbitrary field coefficients. Most applications that use TDA use binary coefficients in the VR filtration, and as a result, most existing TDA software can only handle binary coefficients.  Using other coefficient fields $\mathbb{Z} / p\mathbb{Z}$ for any prime $p$, makes it possible to detect "twists" in point cloud data.  We have included an example notebook that shows how this can differentiate a point cloud sampled from the boundary of a Moebius strip from a point cloud sampled from an ordinary, untwisted loop.

Applications of this surprisingly appears in some important real world scenarios such as periodic time series analysis [@perea2015sliding] and image analysis [@perea2014klein]. This feature of Ripser.py has been used in a pipeline to synthesize slow motion videos [@tralie2018slomoloops] and is currently being used to quantify periodicities in repetitive motions with children with autism spectrum disorder [@tralie2018autism].

## Representative Cocycles
TDA is generally used to quantify topological features, but there has been some research on localizing topological features back in the point cloud. Ripser.py can return "representative cocycles" associated to different homology classes (topological features), for sparse and dense filtrations over any field. Use of this feature can be viewed as topological nonlinear dimension reduction (NLDR). Examples of this can be seen in mapping the point cloud to a circle [@de2011persistent], which is useful for parameterizing periodic data, or in mapping the point cloud to the projective plane [@perea2018multiscale], which shows up in analysis of image patches.

## Higher Order Homology
Because of long standing speed concerns with TDA, most applications have focused on 0-dimensional (connected components) and 1-dimensional (loops) homology features from VR filtrations.  However, Ripser.py is fast enough to compute 2-dimensional homology (voids e.g. empty space in a basketball) for modest point clouds.  This feature has been used to quantify quasiperiodic phenomena in videos of patients with vocal fold disorders [@tralie2017quasi], and we anticipate more researchers will use the feature for other problems now that it is computationally accessible.


# Source Code
The source code for Ripser.py is available on Github through the Scikit-TDA organization [https://github.com/scikit-tda/Ripser.py](https://github.com/scikit-tda/Ripser.py).   The original Ripser library can be found at [https://github.com/Ripser/Ripser/](https://github.com/Ripser/Ripser/)


# Acknowledgements

Christoher Tralie and Rann Bar-On were supported by an NSF big data grant DKA-1447491. Nathaniel Saul was partially supported by NSF DBI-1661348 and by Washington NASA Space Grant Consortium, NASA Grant #NNX15AJ98H. We thank Ulrich Bauer for the original Ripser library and for valuable feedback during development of Ripser.py.  We also thank Jose Perea, William Guss, and Matija Čufar for helpful feedback and bug fixes.  Finally, we thank the students of the "Topological Data Analysis and Persistent Homology" workshop in Levico, Italy for beta testing the code.

# References
