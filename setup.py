from setuptools import setup, Extension
import pybind11
import sys
import os

compile_args = ["/std:c++17", "/O3"] if sys.platform == "win32" else ["-std=c++17", "-O3"]
kinship_core_path = os.path.join("src", "lineagekit", "core", "kinship_core")
source_files = [os.path.join(kinship_core_path, "kinship.cpp")]
library_paths = [kinship_core_path]
include_dirs = [pybind11.get_include()]
include_dirs.extend(library_paths)

ext_modules = [
    Extension(
        "lineagekit.kinship",
        source_files,
        include_dirs=include_dirs,
        language="c++",
        extra_compile_args=compile_args,
    )
]

setup(ext_modules=ext_modules)
