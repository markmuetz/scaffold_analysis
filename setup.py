#!/usr/bin/env python
import os

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

from scaffold.version import get_version


def read(fname):
    try:
        return open(os.path.join(os.path.dirname(__file__), fname)).read()
    except:
        return ''


setup(
    name='scaffold',
    version=get_version(),
    description='High-resolution analysis',
    long_description=read('readme.rst'),
    author='Mark Muetzelfeldt',
    author_email='m.muetzelfeldt@pgr.reading.ac.uk',
    maintainer='Mark Muetzelfeldt',
    maintainer_email='m.muetzelfeldt@pgr.reading.ac.uk',
    packages=[
        'scaffold',
        'scaffold.cycle',
        'scaffold.expt',
        'scaffold.suite',
        'scaffold.tests'
    ],
    scripts=[ ],
    python_requires='>=3.6',
    install_requires=[
        'omnium>=0.10.2',
        'cloud_tracking',
        'f90nml',
        # 'iris',
        'matplotlib',
        'numpy',
        'scipy',
    ],
    package_data={ },
    url='https://github.com/markmuetz/scaffold_analysis',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: C',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        ],
    keywords=[''],
    )
