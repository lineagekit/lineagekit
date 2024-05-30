from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "kinship",
        ["kinship.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["/std:c++17", "/O2"]
    ),
]

setup(
    name="kinship",
    version="0.1",
    ext_modules=ext_modules,
    install_requires=["pybind11"],
)
