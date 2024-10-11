[![DOI](http://joss.theoj.org/papers/10.21105/joss.00925/status.svg)](https://doi.org/10.21105/joss.00925)
[![PyPI version](https://badge.fury.io/py/ripser.svg)](https://badge.fury.io/py/ripser)
[![Downloads](https://img.shields.io/pypi/dm/ripser)](https://pypi.python.org/pypi/ripser/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/ripser.svg)](https://anaconda.org/conda-forge/ripser)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/ripser.svg)](https://anaconda.org/conda-forge/ripser)

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

Ripser.py is available on `pypi` with wheels for all major platforms. To install, type the following command into your environment:

```bash
pip install ripser
```

### Local build

If the above command fails or if you want to develop and contribute to
`ripser.py`, you can build `ripser.py` locally. To do so, clone this
repository. From within the cloned repository, execute `pip install .` to build
locally, or `pip install -e .` for a local,
[editable](https://setuptools.pypa.io/en/latest/userguide/development_mode.html)
build. Either of the above two commands will install all required dependencies.
Explicitly, the dependencies of `ripser.py` are

- Cython,
- numpy,
- scipy,
- scikit-learn,
- persim,

and their required dependencies.

**Windows users:** If you are using a Windows machine, you _may_ also need to install [MinGW](http://www.mingw.org) on your system.

**Mac users:** Updating your Xcode and Xcode command line tools will probably fix any issues you have with installation.

#### Optional dependency

Ripser.py when compiled from source can have a _steroid_<sup>1</sup> shot by replacing the standard `unordered_map` from the STL by one of the fastest implementation available: [robin_hood](https://github.com/martinus/robin-hood-hashing). Benchmarking of Ripser.py using the `robin_hood` implementation showed speed-ups up to **30%**.

To be able to use `robin_hood` instead of STL, you only need to clone the repository containing the implementation:

```
# Run this command at the root of the project
git clone https://github.com/martinus/robin-hood-hashing robinhood
```

After cloning robinhood with the above command, install `ripser.py` with

```
pip install -v .
```

This will install a local version of `ripser.py` with verbose output. In the verbose output,
you will see confirmation that robinhood was found or not.

<sup>1</sup> The Python package is already compiled with `robin_hood` by default.

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

If you use this package, please site the JoSS paper found here [![DOI](http://joss.theoj.org/papers/10.21105/joss.00925/status.svg)](https://doi.org/10.21105/joss.00925) and the JACT paper about Ripser found here [![DOI:10.1007/s41468-021-00071-5](https://zenodo.org/badge/DOI/10.1007/s41468-021-00071-5.svg)](https://doi.org/10.1007/s41468-021-00071-5).

You can use the following bibtex entries:

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

@article{Bauer2021Ripser,
    AUTHOR = {Bauer, Ulrich},
     TITLE = {Ripser: efficient computation of {V}ietoris-{R}ips persistence
              barcodes},
   JOURNAL = {J. Appl. Comput. Topol.},
  FJOURNAL = {Journal of Applied and Computational Topology},
    VOLUME = {5},
      YEAR = {2021},
    NUMBER = {3},
     PAGES = {391--423},
      ISSN = {2367-1726},
   MRCLASS = {55N31 (55-04)},
  MRNUMBER = {4298669},
       DOI = {10.1007/s41468-021-00071-5},
       URL = {https://doi.org/10.1007/s41468-021-00071-5},
}
```
