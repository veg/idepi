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
      packages=[
            'idepi',
            'idepi.alphabet',
            'idepi.argument',
            'idepi.datasource',
            'idepi.filter',
            'idepi.future',
            'idepi.hmmer',
            'idepi.labeler',
            'idepi.logging',
            'idepi.mlpy',
            'idepi.phylogeny',
            'idepi.results',
            'idepi.scripts',
            'idepi.simulation',
            'idepi.svm',
            'idepi.test',
            'idepi.util'
      ],
      package_dir={
            'idepi': 'lib/idepi',
            'idepi.alphabet': 'lib/idepi/alphabet',
            'idepi.argument': 'lib/idepi/argument',
            'idepi.datasource': 'lib/idepi/datasource',
            'idepi.filter': 'lib/idepi/filter',
            'idepi.future': 'lib/idepi/future',
            'idepi.hmmer': 'lib/idepi/hmmer',
            'idepi.labeler': 'lib/idepi/labeler',
            'idepi.logging': 'lib/idepi/logging',
            'idepi.mlpy': 'lib/idepi/mlpy',
            'idepi.phylogeny': 'lib/idepi/phylogeny',
            'idepi.results': 'lib/idepi/results',
            'idepi.scripts': 'lib/idepi/scripts',
            'idepi.simulation': 'lib/idepi/simulation',
            'idepi.svm': 'lib/idepi/svm',
            'idepi.test': 'lib/idepi/test',
            'idepi.util': 'lib/idepi/util',
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
