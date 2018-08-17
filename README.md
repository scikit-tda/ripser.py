[![PyPI version](https://badge.fury.io/py/ripser.svg)](https://badge.fury.io/py/ripser)
[![Build Status](https://travis-ci.org/scikit-tda/ripser.py.svg?branch=master)](https://travis-ci.org/scikit-tda/ripser.py)
[![Build status](https://ci.appveyor.com/api/projects/status/sfy7yybs66e5qanu?svg=true)](https://ci.appveyor.com/project/scikit-tda/ripser)
[![codecov](https://codecov.io/gh/scikit-tda/ripser.py/branch/master/graph/badge.svg)](https://codecov.io/gh/scikit-tda/ripser.py)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

# Ripser

Ripser.py is a renovated Python implementation of the Ripser package. We provide a slim interface for computing persistence cohomology of sparse and dense data sets, visualizing persistence diagrams, computing lowerstar filtrations on images, and computing representative cochains.

To aid your use of the package, we've put together a large set of notebooks that demonstrate many of the features available. Complete documentation about the package can be found at [ripser.scikit-tda.org](https://ripser.scikit-tda.org). 

Through extensive testing and continuous integration, Ripser.py is easy to install on Mac, Linux, and Windows platforms.

If you're looking for the original C++ library, you can find it at [Ripser/ripser](https://github.com/ripser/ripser).

## Setup

Installation requires Cython, but other than that, it is available on all platforms (if you are having trouble installing, please let us know!)

```
pip install Cython
pip install Ripser
```

## Usage

The interface is as simple as can be:

```
import numpy as np
from ripser import ripser, plot_dgms

data = np.random.random((100,2))
diagrams = ripser(data)['dgms']
plot_dgms(diagrams)
```

We also supply a Scikit-learn transformer style object if you would prefer to use that:

```
import numpy as np
from ripser import Rips

rips = Rips()
data = np.random.random((100,2))
diagrams = rips.fit_transform(data)['dgms']
rips.plot(diagrams)
```

# License

Ripser.py is available under an MIT license!

# Contributions

We welcome all kinds of contributions! There are lots of opportunities for potential projects, so please get in touch if you would like to help out. Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.
