"""
Setup script for DNA breathing dynamics framework.

Handles C/C++ extension building with Apple Silicon optimization.
"""

import os
import platform
import sys
from pathlib import Path

from setuptools import Extension, setup

# Detect Apple Silicon
IS_APPLE_SILICON = (
    platform.system() == "Darwin" and platform.machine() == "arm64"
)

# Determine Homebrew prefix
if IS_APPLE_SILICON:
    HOMEBREW_PREFIX = "/opt/homebrew"
else:
    HOMEBREW_PREFIX = "/usr/local"

# Compiler flags
EXTRA_COMPILE_ARGS = [
    "-O3",
    "-march=native",
    "-ffast-math",
    "-Wall",
    "-Wextra",
]

EXTRA_LINK_ARGS = []

# Apple Silicon specific optimizations
if IS_APPLE_SILICON:
    EXTRA_COMPILE_ARGS.extend([
        "-mcpu=apple-m1",  # or detected M2/M3
        "-DAPPLE_SILICON",
    ])

# Include and library paths
INCLUDE_DIRS = [
    "src/extensions",
    f"{HOMEBREW_PREFIX}/include",
]

LIBRARY_DIRS = [
    f"{HOMEBREW_PREFIX}/lib",
]

LIBRARIES = [
    "mpfr",
    "gmp",
]

# Define extensions (currently empty, ready for future C/C++ modules)
extensions = [
    # Example extension (uncomment when needed):
    # Extension(
    #     "src.extensions.fast_compute",
    #     sources=["src/extensions/fast_compute.c"],
    #     include_dirs=INCLUDE_DIRS,
    #     library_dirs=LIBRARY_DIRS,
    #     libraries=LIBRARIES,
    #     extra_compile_args=EXTRA_COMPILE_ARGS,
    #     extra_link_args=EXTRA_LINK_ARGS,
    # ),
]

if __name__ == "__main__":
    setup(
        ext_modules=extensions,
    )
