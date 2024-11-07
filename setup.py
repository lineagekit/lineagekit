import os
import sys
from setuptools import setup, Extension, find_packages
import pybind11


def parse_requirements(filename):
    with open(filename, 'r') as f:
        return f.read().splitlines()


# Determine the appropriate compiler flags based on the operating system
if sys.platform == "win32":
    compile_args = ["/std:c++17", "/O3"]
else:
    compile_args = ["-std=c++17", "-O3"]


# Define the extension modules
def get_extensions():
    ext_modules = []
    kinship_core_dir = os.path.join('src', 'basic', 'kinship_core')

    ext_modules.append(
        Extension(
            "kinship",
            [os.path.join(kinship_core_dir, "kinship.cpp")],  # C++ source file
            include_dirs=[pybind11.get_include()],  # Include directories for pybind11
            language="c++",
            extra_compile_args=compile_args,  # Compiler flags
        )
    )

    return ext_modules


# Setup configuration
setup(
    name="lineagekit",
    version="1.0",
    packages=find_packages(where='src'),  # Find all packages under 'src'
    package_dir={'': 'src'},  # Set the package root directory
    ext_modules=get_extensions(),
    install_requires=parse_requirements('requirements.txt'),
    package_data={'': ['kinship.pyi']},
    zip_safe=False,
)
