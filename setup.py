#!/usr/bin/env python

from __future__ import division, print_function

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from idepi import __version__ as _idepi_version

setup(name='idepi',
      version=_idepi_version,
      description='IDentify EPItope',
      author='N Lance Hepler',
      author_email='nlhepler@gmail.com',
      url='http://github.com/veg/idepi',
      license='GNU GPL version 3',
      packages=['idepi', 'idepi._mlpy', 'idepi._scripts'],
      package_dir={
            'idepi': 'lib/idepi',
            'idepi._mlpy': 'lib/idepi/_mlpy',
            'idepi._scripts': 'lib/idepi/_scripts'
      },
      package_data={'idepi': ['data/*', 'hyphy/*.bf']},
      scripts=['scripts/idepi'],
      requires=[
            'Bio (>=1.58)',
            'BioExt (>=0.9.4)',
            'cvxopt (>=1.1.5)',
            'fakemp (>=0.9.1)',
            'hppy (>=0.9.2)',
            'mlpy (>=3.5.0)',
            'mrmr (>=0.9.2)',
            'numpy (>=1.6)',
            'pyxval (>=0.9.3)',
            'six'
      ]
     )
