.. ripser documentation master file, created by
   sphinx-quickstart on Sun Jul 22 20:37:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|PyPI version| |Travis-CI| |Appveyor| |Codecov| |License: MIT|

Ripser.py
==========

Ripser.py is a renovated Python implementation of the Ripser package. 

We provide a slim interface for computing persistence cohomology of sparse and dense data sets, visualizing persistence diagrams, computing lowerstar filtrations on images, and computing representative cochains. 

We supply a large set of interactive notebooks that demonstrate how to take advantage of all the features available.

Through extensive testing and continuous integration, Ripser.py is easy to install on Mac, Linux, and Windows platforms.   

You can find the source code on github at `Scikit-TDA/ripser.py <https://github.com/scikit-tda/ripser.py>`_. For the original C++ library, see `Ripser/ripser <https://github.com/Ripser/ripser/releases/latest>`_.

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
    Basic Usage

.. toctree::
    :maxdepth: 2
    :caption: Tutorials

    Representative Cocycles
    Approximate Sparse Filtrations
    parse Distance Matrices
    Lower Star Time Series
    Lower Star Image Filtrations
    Moebius Strip And The Field of Coefficients

.. toctree::
    :maxdepth: 2
    :caption: API Reference
    
    reference


.. |PyPI version| image:: https://badge.fury.io/py/ripser.svg
   :target: https://badge.fury.io/py/ripser

.. |Travis-CI| image:: https://travis-ci.org/scikit-tda/ripser.py.svg?branch=master
    :target: https://travis-ci.org/scikit-tda/ripser.py

.. |Appveyor| image:: https://ci.appveyor.com/api/projects/status/sfy7yybs66e5qanu?svg=true
    :target: https://ci.appveyor.com/project/scikit-tda/ripser.py
.. |Codecov| image:: https://codecov.io/gh/scikit-tda/ripser.py/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/scikit-tda/ripser.py
.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT)
