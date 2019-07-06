#!/usr/bin/env python3

import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name="graph-tools",
    version="1.0",
    author="Hiroyuki Ohsaki",
    author_email="ohsaki@lsnl.jp",
    description=
    "tools for graph theory and network science with many generation models",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/h-ohsaki/graph-tools",
    packages=setuptools.find_packages(),
    install_requires=['numpy', 'perlcompat', 'pytess'],
#    scripts=[''],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
)
