---
title: 'Ripser.py: A Lean Persistent Homology Library for Python'
tags:
  - Python
  - Persistent Homology
  - Topological Data Analysis
authors:
  - name: Christopher Tralie
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Nathaniel Saul
    orcid: 0000-0000-0000-0000
    affiliation: 2
affiliations:
 - name: Department of Computer Science, Duke University
   index: 1
 - name: Department of Mathematics and Statistics, Washington State University
   index: 2
date: 15 August 2018
bibliography: paper.bib
---

# Notes

- This paper was put together by following directions from https://joss.readthedocs.io/en/latest/submitting.html

# Summary

We provide a renovated Python implementation of the Ripser package. 

Using Cython, we have rebuilt the Ripser C++ code to provide an optimized Python interface. 

Because of the unprecidented speed of Ripser, it has created a large user base for both research and applications. The library as stands is only as a command line tool. Multiple efforts have been made to wrap the C++ library for use in other languages.  Our implementation is unique because, by working closely with the original Ripser developer Ulrich Baeur, Ripser.py provides access to cocycles for each of the cohomological generators, . 

This Python library supplies two interfaces, one lightweight and functional interface as well as an object-oriented interface designed to fit within the Scikit-Learn transformer paradigm [@scikit-learn]. 

Unique to Ripser.py is functionality for computing sublevelset filtrations on images. Sublevelset filtrations have been shown to be strongly applicable to descriminating images with nuanced differences [@obayashi2018persistence]

The source code for Ripser.py is available on github through the Scikit-TDA organization [https://github.com/scikit-tda/ripser.py][https://github.com/scikit-tda/ripser.py]. 


# Acknowledgements

We acknowledge contributions from [Everyone who contributed to the GH repo, Ulrich Baeur, others?]

# References