.. ripser documentation master file, created by
   sphinx-quickstart on Sun Jul 22 20:37:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|PyPI version| |Travis-CI| |Appveyor| |Codecov| |License: MIT|

Ripser.py
==========

Ripser.py 

This package provides the awesome Ripser project as an easy to use Python module. It is easy to install and is even easier to use.

For the original C++ library, see `Ripser/ripser <https://github.com/Ripser/ripser/releases/latest>`_.

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
    :target: https://ci.appveyor.com/project/scikit-tda/ripser
.. |Codecov| image:: https://codecov.io/gh/scikit-tda/ripser.py/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/scikit-tda/ripser.py
.. |License: MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT)
