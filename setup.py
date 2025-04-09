#!/bin/envs/python3

from setuptools import setup, Extension
import numpy as np

if __name__ == "__main__":
    np_inc = np.get_include()
    setup(
        name="ARTmie",
        version="0.1.0",
        description="Fast Mie calculation library with C++ backend",
        authors=["Enrico P. Metzner"],
        ext_modules=[
            Extension(
                'ARTmie',
                sources=['src/ARTmie.cpp'],
                include_dirs=[np_inc]
            ),
        ]
    )
