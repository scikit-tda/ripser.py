import sys
import os
import platform

from setuptools import setup
from setuptools.extension import Extension

# Ensure Cython is installed before we even attempt to install Ripser.py
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org or install it with `pip install Cython`")
    sys.exit(1)

## Get version information from _version.py
import re
VERSIONFILE="ripser/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

# Use README.md as the package long description  
with open('README.md') as f:
    long_description = f.read()

class CustomBuildExtCommand(build_ext):
    """ This extension command lets us not require numpy be installed before running pip install ripser 
        build_ext command for use when numpy headers are needed.
    """

    def run(self):
        # Import numpy here, only when headers are needed
        import numpy
        # Add numpy headers to include_dirs
        self.include_dirs.append(numpy.get_include())
        # Call original build_ext command
        build_ext.run(self)

extra_compile_args = ["-Ofast", "-D_hypot=hypot"]
extra_link_args = []

if platform.system() == "Windows":
    extra_compile_args.extend([
        # Supported by Visual C++ >=14.1
        '/std:c++14'
    ])
elif platform.system() == "Darwin":
    extra_compile_args.extend([
        '-std=c++11', 
        "-mmacosx-version-min=10.9"
    ])
    extra_link_args.extend([
        "-stdlib=libc++", 
        "-mmacosx-version-min=10.9"
    ])
else:
    extra_compile_args.extend([
        "-std=c++11"
    ])

macros = [
    ("USE_COEFFICIENTS", 1),
    ("NDEBUG", 1),
    ("ASSEMBLE_REDUCTION_MATRIX", 1)
]

# Robinhood
robinhood_path = os.path.join('ripser', 'robinhood')
if os.path.isdir(robinhood_path):
    macros.extend([("USE_ROBINHOOD_HASHMAP", 1)])

    robinhood_include_path = os.path.join('src', 'include')
    if platform.system() == "Windows":
        extra_compile_args.extend(['/I'+os.path.join(robinhood_path,
                                                     robinhood_include_path)])
    else:
        extra_compile_args.extend(['-I'+os.path.join(robinhood_path,
                                                     robinhood_include_path)])

ext_modules = Extension(
    "pyRipser",
    sources=["ripser/pyRipser.pyx"],
    define_macros=macros,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    language="c++"
)


setup(
    name="ripser",
    version=verstr,
    description="A Lean Persistent Homology Library for Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Chris Tralie, Nathaniel Saul",
    author_email="chris.tralie@gmail.com, nat@riverasaul.com",
    url="https://ripser.scikit-tda.org",
    license='MIT',
    packages=['ripser'],
    ext_modules=cythonize(ext_modules),
    install_requires=[
        'Cython',
        'numpy',
        'scipy',
        'scikit-learn',
        'persim'
    ],
    extras_require={
        'testing': [ # `pip install -e ".[testing]"``
            'pytest'  
        ],
        'docs': [ # `pip install -e ".[docs]"`
            'sktda_docs_config'
        ],
        'examples': [
            'persim',
            'tadasets',
            'jupyter',
            'pillow'
        ]
    },
    cmdclass={'build_ext': CustomBuildExtCommand},
    python_requires='>=3.6',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Financial and Insurance Industry',
        'Intended Audience :: Healthcare Industry',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8'
    ],
    keywords='persistent homology, rips filtration, persistence diagrams, topology data analysis, algebraic topology, unsupervised learning'
)
