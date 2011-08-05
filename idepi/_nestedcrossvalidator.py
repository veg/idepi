#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from copy import deepcopy

from _crossvalidator import CrossValidator
from _gridsearcher import GridSearcher
from _perfstats import PerfStats


__all__ = ['NestedCrossValidator']


class NestedCrossValidator(CrossValidator):

    ACCURACY            = PerfStats.ACCURACY
    PPV, PRECISION      = PerfStats.PPV, PerfStats.PRECISION
    NPV                 = PerfStats.NPV
    SENSITIVITY, RECALL = PerfStats.SENSITIVITY, PerfStats.RECALL
    SPECIFICITY, TNR    = PerfStats.SPECIFICITY, PerfStats.TNR
    FSCORE              = PerfStats.FSCORE
    MINSTAT             = PerfStats.MINSTAT

    def __init__(self, classifiercls, folds, cv={}, mode=None, optstat=PerfStats.MINSTAT, gs={}):

        ncv = {
            'classifiercls': classifiercls,
            'folds': folds - 1,
            'optstat': optstat,
            'cv': cv,
            'mode': mode,
            'gs': gs,
        }

        super(NestedCrossValidator, self).__init__(GridSearcher, folds, ncv, mode)
