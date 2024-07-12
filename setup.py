import os
import platform
import re
import sys

from setuptools import Extension, setup

# Ensure Cython is installed before we even attempt to install Ripser.py
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org or install it with `pip install Cython`")
    sys.exit(1)


def get_version():
    VERSIONFILE = "src/ripser/_version.py"
    verstrline = open(VERSIONFILE, "rt").read()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        return mo.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


class CustomBuildExtCommand(build_ext):
    """This extension command lets us not require numpy be installed before running pip install ripser
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
    extra_compile_args.extend(
        [
            # Supported by Visual C++ >=14.1
            "/std:c++14"
        ]
    )
elif platform.system() == "Darwin":
    extra_compile_args.extend(["-std=c++11", "-mmacosx-version-min=10.9"])
    extra_link_args.extend(["-stdlib=libc++", "-mmacosx-version-min=10.9"])
else:
    extra_compile_args.extend(["-std=c++11"])

macros = [("USE_COEFFICIENTS", 1), ("NDEBUG", 1), ("ASSEMBLE_REDUCTION_MATRIX", 1)]

# Robinhood
robinhood_path = os.path.join("robinhood")
if os.path.isdir(robinhood_path):
    print(
        "\nFound local copy of robinhood! Using robinhood for ripser.py compilation.\n"
    )
    macros.extend([("USE_ROBINHOOD_HASHMAP", 1)])

    robinhood_include_path = os.path.join("src", "include")
    if platform.system() == "Windows":
        extra_compile_args.extend(
            ["/I" + os.path.join(robinhood_path, robinhood_include_path)]
        )
    else:
        extra_compile_args.extend(
            ["-I" + os.path.join(robinhood_path, robinhood_include_path)]
        )
else:
    print("Did not find a local copy of robinhood. Proceeding anyways.")

ext_modules = Extension(
    "pyRipser",
    sources=["src/ripser/pyRipser.pyx"],
    define_macros=macros,
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    language="c++",
)


setup(
    ext_modules=cythonize(ext_modules),
    cmdclass={"build_ext": CustomBuildExtCommand},
    version=get_version(),
)
