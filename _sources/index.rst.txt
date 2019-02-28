.. ripser documentation master file, created by
   sphinx-quickstart on Sun Jul 22 20:37:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|DOI| |PyPI version| |Downloads| |Travis-CI| |Appveyor| |Codecov| |License: MIT|

Ripser.py
==========

Ripser.py is a lean persistent homology package for Python. Building on the blazing fast C++ Ripser package as the core computational engine, Ripser.py provides an intuitive interface for 

- computing persistence cohomology of sparse and dense data sets, 
- visualizing persistence diagrams, 
- computing lowerstar filtrations on images, and 
- computing representative cochains. 


Additionally, through extensive testing and continuous integration, Ripser.py is easy to install on Mac, Linux, and Windows platforms.

We supply a large set of interactive notebooks that demonstrate how to take advantage of all the features available.

You can find the source code on github at `Scikit-TDA/Ripser.py <https://github.com/scikit-tda/ripser.py>`_. For the original C++ library, see `Ripser/ripser <https://github.com/Ripser/ripser/releases/latest>`_.


Example Usage
--------------



.. code:: python

    import numpy as np
    from ripser import Rips
    r = Rips()

    data = np.random.random((100,2))
    diagram = r.fit_transform(data)
    r.plot(diagram, show=True)

.. toctree::
    :maxdepth: 2
    :caption: Background

    about
    notebooks/Basic Usage
    reference/index

.. toctree::
    :maxdepth: 2
    :caption: Tutorials

    notebooks/Representative Cocycles
    notebooks/Approximate Sparse Filtrations
    notebooks/Sparse Distance Matrices
    notebooks/Lower Star Time Series
    notebooks/Lower Star Image Filtrations
    notebooks/Moebius Strip And The Field of Coefficients
    notebooks/Greedy Subsampling for Fast Approximate Computation




.. |DOI| image:: http://joss.theoj.org/papers/10.21105/joss.00925/status.svg
    :target: https://doi.org/10.21105/joss.00925

.. |Downloads| image:: https://pypip.in/download/ripser/badge.svg
    :target: https://pypi.python.org/pypi/ripser/

.. |PyPI version| image:: https://badge.fury.io/py/ripser.svg
   :target: https://badge.fury.io/py/ripser

.. |Travis-CI| image:: https://travis-ci.org/scikit-tda/ripser.py.svg?branch=master
    :target: https://travis-ci.org/scikit-tda/ripser.py

.. |Appveyor| image:: https://ci.appveyor.com/api/projects/status/020nrvrq2rdg2iu1?svg=true
    :target: https://ci.appveyor.com/project/sauln/ripser-py

.. |Codecov| image:: https://codecov.io/gh/scikit-tda/ripser.py/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/scikit-tda/ripser.py

.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT)
