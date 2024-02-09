"""A setuptools based setup module.

See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name="maskGenerator",  # Required

    version="0.1.0",  # Required

    description="functions and scripts for masking and normalising BM31 total scattering data",  # Optional

    #url="https://github.com/pypa/sampleproject",  # Optional
    author="K. P. Marshall",  # Optional
    # This should be a valid email address corresponding to the author listed
    # above.
    author_email="kenneth.marshall@esrf.fr",  # Optional

    # When your source code is in a subdirectory under the project root, e.g.
    # `src/`, it is necessary to specify the `package_dir` argument.
    package_dir={'maskGenerator':'.'},
    packages=['maskGenerator'],  # Required
    install_requires=['fabio','pyfai', 'cryio','matplotlib'],
    python_requires=">=3.7",
)
