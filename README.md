[![DOI](http://joss.theoj.org/papers/10.21105/joss.00925/status.svg)](https://doi.org/10.21105/joss.00925)
[![PyPI version](https://badge.fury.io/py/ripser.svg)](https://badge.fury.io/py/ripser)
[![Downloads](https://pypip.in/download/ripser/badge.svg)](https://pypi.python.org/pypi/ripser/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/ripser.svg)](https://anaconda.org/conda-forge/ripser)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/ripser.svg)](https://anaconda.org/conda-forge/ripser)

[![Build Status](https://travis-ci.org/scikit-tda/ripser.py.svg?branch=master)](https://travis-ci.org/scikit-tda/ripser.py)
[![Build status](https://ci.appveyor.com/api/projects/status/020nrvrq2rdg2iu1?svg=true)](https://ci.appveyor.com/project/sauln/ripser-py)
[![codecov](https://codecov.io/gh/scikit-tda/ripser.py/branch/master/graph/badge.svg)](https://codecov.io/gh/scikit-tda/ripser.py)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# Ripser.py

Ripser.py is a lean persistent homology package for Python. Building on the blazing fast C++ Ripser package as the core computational engine, Ripser.py provides an intuitive interface for 

- computing persistence cohomology of sparse and dense data sets, 
- visualizing persistence diagrams, 
- computing lowerstar filtrations on images, and 
- computing representative cochains. 

Additionally, through extensive testing and continuous integration, Ripser.py is easy to install on Mac, Linux, and Windows platforms.

To aid your use of the package, we've put together a large set of notebooks that demonstrate many of the features available. Complete documentation about the package can be found at [ripser.scikit-tda.org](https://ripser.scikit-tda.org). 


## Related Projects

If you're looking for the original C++ library, you can find it at [Ripser/ripser](https://github.com/ripser/ripser).

If you're looking for a GPU-accelerated version of Ripser, you can find it at [Ripser++](https://github.com/simonzhang00/ripser-plusplus)

## Setup


Ripser.py is available on all major platforms. All that is required is that you install the standard Python numerical computing libraries and Cython. 

Dependencies:
- Cython
- numpy
- scipy
- scikit-learn
- persim

**Windows users:** If you are using a Windows machine, you will also need to install [MinGW](http://www.mingw.org) on your system.

**Mac users:** Updating your Xcode and Xcode command line tools will probably fix any issues you have with installation.

Cython should be the only library required before installation.  To install, type the following commands into your environment:

```
pip install cython
pip install ripser
```

If you are having trouble installing, please let us know!


## Usage

The interface is as simple as can be:

```
import numpy as np
from ripser import ripser
from persim import plot_diagrams

data = np.random.random((100,2))
diagrams = ripser(data)['dgms']
plot_diagrams(diagrams, show=True)
```

We also supply a Scikit-learn transformer style object if you would prefer to use that:

```
import numpy as np
from ripser import Rips

rips = Rips()
data = np.random.random((100,2))
diagrams = rips.fit_transform(data)
rips.plot(diagrams)
```

<img src="https://i.imgur.com/WmQPYnn.png" alt="Ripser.py output persistence diagram" width="70%"/>

# Contributions

We welcome all kinds of contributions! Please get in touch if you would like to help out. Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

If you found a bug, have questions, or are just having trouble with the library, please open an issue in our [issue tracker](https://github.com/scikit-tda/ripser.py/issues/new) and we'll try to help resolve the concern.

# License

Ripser.py is available under an MIT license! The core C++ code is derived from Ripser, which is also available under an MIT license and copyright to Ulrich Bauer. The modifications, Python code, and documentation is copyright to Christopher Tralie and Nathaniel Saul.

# Citing

If you use this package, please site the JoSS paper found here: [![DOI](http://joss.theoj.org/papers/10.21105/joss.00925/status.svg)](https://doi.org/10.21105/joss.00925)

You can use the following bibtex entry
```
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
```
