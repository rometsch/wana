#!/usr/bin/env python3
import os
from setuptools import setup, find_namespace_packages

setup(  name="wana"
        ,version="0.1"
        ,description="Walking analysis utility to read and analyze sensor data."
        ,author="Thomas Rometsch"
        ,package_dir={'': 'src'}
        ,packages=find_namespace_packages(where="src")
        ,install_requires=[
          "numpy",
          "matplotlib"
        ]
        ,zip_safe=False
        ,entry_points = {
            'console_scripts': ['wana=wana._command_line_:main'],
        }
)
