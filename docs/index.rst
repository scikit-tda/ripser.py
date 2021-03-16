.. ripser documentation master file, created by
   sphinx-quickstart on Sun Jul 22 20:37:23 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|DOI| |PyPI version| |Downloads| |License: MIT|

|Travis-CI| |Appveyor| |Codecov| 

Ripser.py is a lean persistent homology package for Python. Building on the blazing fast C++ Ripser package as the core computational engine, Ripser.py provides an intuitive interface for 

- computing persistence cohomology of sparse and dense data sets, 
- visualizing persistence diagrams, 
- computing lowerstar filtrations on images, and 
- computing representative cochains. 

We supply a large set of interactive notebooks that demonstrate how to take advantage of all the features available.

Ripser.py is an evolution of the original C++ Ripser project. We have put extensive work into making the package available to Python developers across all major platforms. If you are having trouble installing, please let us know by opening a github issue.


You can find the source code on github at `Scikit-TDA/Ripser.py <https://github.com/scikit-tda/ripser.py>`_. For the original C++ library, see `Ripser/ripser <https://github.com/Ripser/ripser>`_.



Setup
------

Ripser.py is available on Pypi. To install, you'll first need Cython. 

.. code:: python

    pip install Cython
    pip install Ripser

Example Usage
--------------

The interface is as simple as can be:


.. code:: python

    import numpy as np
    from ripser import ripser
    from persim import plot_diagrams

    data = np.random.random((100,2))
    diagrams = ripser(data)['dgms']
    plot_diagrams(diagrams, show=True)


We also supply a Scikit-learn transformer style object if you would prefer to use that:

.. code:: python

    import numpy as np
    from ripser import Rips

    rips = Rips()
    data = np.random.random((100,2))
    diagrams = rips.fit_transform(data)
    rips.plot(diagrams)


.. image:: https://i.imgur.com/WmQPYnn.png


Source Compilation
---------------------

In order to compile, you'll first need to do:

.. code:: bash

    git clone https://github.com/scikit-tda/ripser.py
    cd ripser.py
    pip install -e .

You'll then be able to use Ripser.py as in the previous examples but compiled directly from sources.
In order to obtain the best experience in performances when compiling from sources, you'll need an additional library, robin_hood_hashmap. To be able to use compile with this library enable, do:

.. code:: bash

    git clone https://github.com/martinus/robin-hood-hashing ripser/robinhood

The following table shows a comparison of performances with and without robin_hood_hashmap:

+------------+--------+-------------+-------+---------+--------------+------------------+
| Dataset    | size   | threshold   | dim   | coeff   | normal [s]   | robin_hood [s]   |
+============+========+=============+=======+=========+==============+==================+
| sphere3    | 192    |             | 2     | 2       | 1.5          | 1.2              |
+------------+--------+-------------+-------+---------+--------------+------------------+
| dragon     | 2000   |             | 1     | 2       | 2.9          | 2.5              |
+------------+--------+-------------+-------+---------+--------------+------------------+
| o3         | 1024   | 1.8         | 3     | 2       | 2.9          | 2.2              |
+------------+--------+-------------+-------+---------+--------------+------------------+
| random16   | 50     |             | 7     | 2       | 8.4          | 6.0              |
+------------+--------+-------------+-------+---------+--------------+------------------+
| fractal    | 512    |             | 2     | 2       | 17.7         | 14               |
+------------+--------+-------------+-------+---------+--------------+------------------+
| o3         | 4096   | 1.4         | 3     | 2       | 68.6         | 53.4             |
+------------+--------+-------------+-------+---------+--------------+------------------+

Contributions
----------------

We welcome contributions of all shapes and sizes. There are lots of opportunities for potential projects, so please get in touch if you would like to help out. Everything from an implementation of your favorite distance, notebooks, examples, and documentation are all equally valuable so please don’t feel you can’t contribute.

To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

Citing
----------

If you use this package, please site the JoSS paper found here: |DOI|

You can use the following bibtex entry

::

    @article{ctralie2018ripser,
        doi = {10.21105/joss.00925},
        url = {https://doi.org/10.21105/joss.00925},
        year  = {2018},
        month = {Sep},
        publisher = {The Open Journal},
        volume = {3},
        number = {29},
        pages = {925},
        author = {Christopher Tralie and Nathaniel Saul and Rann Bar-On},
        title = {{Ripser.py}: A Lean Persistent Homology Library for Python},
        journal = {The Journal of Open Source Software}
    }




License
--------

Ripser.py is available under an MIT license! The core C++ code is derived from Ripser, which is also available under an MIT license and copyright to Ulrich Baeur. The modifications, Python code, and documentation is copyright to Christopher Tralie and Nathaniel Saul.




.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Background

    notebooks/Basic Usage
    reference/index

.. toctree::
    :maxdepth: 2
    :hidden:
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
