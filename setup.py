import os
import sys

from setuptools import setup, Extension, find_packages
import pybind11

# Determine the appropriate compiler flags based on the operating system
if sys.platform == "win32":
    compile_args = ["/std:c++17", "/O2"]
else:
    compile_args = ["-std=c++17", "-O2"]


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
    name="kinship",
    version="0.1",
    packages=find_packages(include=['kinship', 'kinship.*']),  # Find all packages in the kinship directory
    ext_modules=get_extensions(),
    install_requires=["pybind11"],  # Ensure pybind11 is installed
    zip_safe=False,
)
