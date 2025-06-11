#!/bin/envs/python3

from setuptools import setup, Extension
import numpy as np
import os

if __name__ == "__main__":
    cwd = os.getcwd()
    print('pwd   ',cwd)
    np_inc = np.get_include()
    print('numpy ',np_inc)
    setup(
        name="ARTmie",
        version="0.1.1",
        description="Fast Mie calculation library with C++ backend",
        authors=["Enrico P. Metzner"],
        ext_modules=[
            Extension(
                'ARTmie',
                sources=['src/ARTmie.cpp'],
                include_dirs=[np_inc,cwd+'/src']
            ),
        ]
    )
