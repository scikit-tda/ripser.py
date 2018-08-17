About
======


Ripser.py is an evolution of the original C++ Ripser project. We have put extensive work into making the package available to Python developers across all major platforms. If you are having trouble installing, please let us know by opening a github issue.


Setup
------

Ripser.py is available on Pypi. To install, you'll first need Cython. 

.. code:: python

    pip install Cython
    pip install Ripser

Usage
------

.. code:: python

    import numpy as np
    from ripser import ripser, plot_dgms

    data = np.random.random((100,2))
    diagrams = ripser(data)['dgms']
    plot_dgms(diagrams, show=True)



Note that there is also a *Rips* object with the same functionality, which conforms to the Scikit-learn API.