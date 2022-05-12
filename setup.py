#!/usr/bin/env python
#-*- coding:utf-8 -*-


from setuptools import setup, find_packages

setup(
    name="sinfonia",
    version="0.0.3",
    keywords=("pip", "sinfonia"),
    description="SINFONIA: scalable identification of spatially variable genes for deciphering spatial domains",
    long_description="SINFONIA provides an effective and efficient way to identify spatially variable genes for deciphering spatial domains. We provide documentation in the form of functional application programming interface documentation, tutorials and example workflows at https://sinfonia-svg.readthedocs.io/en/latest/index.html. All SINFONIA wheels distributed on PyPI are MIT licensed.",
    license="MIT Licence",
    url="https://github.com/BioX-NKU/SINFONIA",
    author="Shengquan Chen",
    packages=find_packages(),
    python_requires='>3.6.0',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[
        'numpy>=1.21.6,<1.22',
        'pandas>=1.4.2',
        'scipy>=1.8.0',
        'scikit-learn>=1.0.2',
        'numba>=0.55.1',
        'scanpy>=1.9.1'
    ]
)