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
  orcid: 0000-0000-0000-0000
  affiliation: 2
affiliations:
 - name: Department of Mathematics, Duke University
  index: 1
 - name: Department of Mathematics and Statistics, Washington State University
  index: 2
date: 15 August 2018
bibliography: paper.bib
---


<!---
- This paper was put together by following directions from https://joss.readthedocs.io/en/latest/submitting.html
TODO:
*Citations Other existing software packages, including julia Ripser wrapper?
-->




# Summary

## Background
Topological data analysis (TDA) [@edelsbrunner2010computational],[@carlsson2009topology] is a field focused on computing topological descriptors of point cloud data, summarizing such features as connected components, loops, and voids, at multiple scales in a data structure known as a "persistence diagram." TDA has found wide applications across nonlinear time series analysis [@perea2015sliding], computer vision [@perea2014klein], computational neuroscience [@giusti2015clique], computational biology [@iyer2010imaging],[@wu2017optimal], and materials science [@kramar2013persistence], to name a few of the many areas of impact in recent years.

The persistence algorithm [@edelsbrunner2010computational] is the main workhorse of TDA, and the most commonly used scheme for generating persistence diagrams is the so-called "Vietoris Rips filtration," since it is easily defined for any point cloud described by pairwise similarites between points. In its naive implementation, it is prohibitively slow, but recently a C++ library known as ``Ripser'' [@bauer2017Ripser] has been devised to aggregate all known computational speedups of the Vietoris Rips filtration into one concise implementation. Because of the unprecidented speed of Ripser, it has created a large user base for both research and applications. However, the library as stands is only as a command line tool. Multiple efforts have been made to wrap the C++ library for use in other languages, often via clunky system calls to the command line from the respective languages. In this work, we provide a renovated Python implementation of the Ripser package using Cython to wrap around the C++ codebase. We have gone through extensive testing via continuous integration frameworks Travis and AppVeyor to ensure it works cross-platform. As a result, it is currently as easy to setup as "pip install Ripser," which helpful for everyone who is interested in TDA, especially for mathematicians with little programming experience who would like to use TDA as an entry point into data science.  We have also created multiple Jupyter notebooks to showcase our API and to hint at some of the many applications that are possible with this library.

## Library Details

Our Python library supplies two interfaces: one lightweight and functional interface, as well as an object-oriented interface designed to fit within the Scikit-Learn transformer paradigm [@scikit-learn]. We have merged together multiple branches of the original Ripser library ("sparse," "representative-cocycles") to expose some lesser known but incredibly useful features hidden in Ripser.  Below we detail some of the special features made easy with our library.

### Sparse filtrations
We allow the user to pass on a sparse distance matrix between points using the scipy.sparse library. In one of our notebooks, we demonstrate Ripser can be used in conjunction with "sparse filtration" approximation algorithms [@cavannageometric] to add further computational speedups.  Additionally, sparse matrices allow us to easily define a "sublevelset filtration," or the "watershed algorithm," for quantifying critical points in grid data such as time series and images, and we have exposed a function to make this easy for image data.  One of our notebooks demonstrates how this can be used to identify cells in an image, for instance.

### Coefficient Fields
Most applications that use TDA use the field $\mathbb{Z} / 2\mathbb{Z}$, or "binary coefficients," in the Rips filtration, and as a result, most TDA software to date can only handle binary coefficients.  However, using other coefficient fields $\mathbb{Z} / p\mathbb{Z}$ for prime $p$, makes it possible to detect "twists" in point cloud data.  This surprisingly appears in some important world applications in periodic time series analysis [@perea2015sliding] and image analysis [@perea2014klein].  Ripser naturally allows the specification of arbitrary field coefficients.  We have included an example notebook that shows how this can differentiate a point cloud sampled from the boundary of a Moebius strip from a point cloud sampled from an ordinary, untwisted loop.  This feature of ripser.py is currently being used to help quantify periodicites in repetitive motions with children with autism spectrum disorder [@tralie2018autism].


### Higher Homology
Because of long standing speed issues with TDA, most applications have focused on 0D (connected components) and 1D (loops) homology features from Rips filtrations.  However, Ripser is fast enough to compute 2D homology for modest point clouds, which enables detection and quantification of voids, so we naturally support this in Ripser.py.  This feature has been used to quantify quasiperiodic phenomena in videos of patients with vocal fold disorders [@tralie2017quasi], and we anticipate more researchers will use the feature for other problems now that it is computationally accessible.

### Representative Cocycles
TDA is generally used to quantify topological features, but there has been some research on localizing topological features back in the point cloud, which can be viewed as topological nonlinear dimension reduction (NLDR).  Ripser.py can return "representative cocycles" associated to different homology classes (topological features), for sparse and dense filtrations over any field.  This is an important step in both circular coordinates, or mapping the point cloud to a circle [@de2011persistent], which is useful for parameterizing periodic data.  It also enables an implementation of projective coordinates [@perea2018multiscale], or mapping the point cloud to the projective plane, which shows up in analysis of image patches.  To date, there has only been one end to end implementation of circular coordinates and zero implementations of projective coordinates, but Ripser.py has enabled the development of both in an offshoot project called "DREiMac" (Dimension Reduction with Eilenberg-MacClane coordinates), which we hope can be used as blackbox tools for understanding circular and twisted features of point clouds in high dimensions.


### Source Code
The source code for Ripser.py is available on github through the Scikit-TDA organization [https://github.com/scikit-tda/Ripser.py][https://github.com/scikit-tda/Ripser.py].   The original Ripser library can be found at [https://github.com/Ripser/Ripser/][https://github.com/Ripser/Ripser/]


# Acknowledgements

We thank Ulrich Bauer for the original Ripser library and for valuable feedback during development of Ripser.py.  We also thank various   Finally, we thank the students of the "Topological Data Analysis and Persistent Homology" workshop in Levico, Italy for beta testing the code.

# References
