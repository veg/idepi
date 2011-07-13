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
from itertools import chain
from math import floor
from random import shuffle
from types import FunctionType

import numpy as np

from _perfstats import PerfStats
from _util import ystoconfusionmatrix


__all__ = ['CrossValidator']


# implement cross-validation interface here, grid-search optional
class CrossValidator(object):

    def __init__(self, classifiercls, folds, cv={}):
        self.classifiercls = classifiercls
        self.folds = folds
        self.cv = cv
        classifiercls_dir = dir(classifiercls)
        self.__learnfunc = None
        for m in ('learn', 'compute'):
            if m in classifiercls_dir:
                self.__learnfunc = m
                break
        if self.__learnfunc is None:
            raise ValueError('No known learning mechanism in base class `%s\'' % repr(classifiercls))
        self.__predictfunc = None
        for m in ('predict', 'pred'):
            if m in classifiercls_dir:
                self.__predictfunc = m
                break
        if self.__predictfunc is None:
            raise ValueError('No known prediction mechanism in base class `%s\'' % repr(classifiercls))

    @staticmethod
    def __partition(l, folds):
        npf = int(floor(l / folds)) # num per fold
        r = l % folds
        p = list(chain(*[[i] * npf for i in xrange(folds)])) + range(r)
        shuffle(p)
        assert(len(p) == l)
        return p

    def crossvalidate(self, x, y, cv={}, extra=None):
        if extra is not None:
            if not isinstance(extra, str) and not isinstance(extra, FunctionType):
                raise ValueError('the `extra\' argument takes either a string or a function.')

        kwargs = deepcopy(self.cv)
        kwargs.update(cv)

        nrow, _ = x.shape

        p = CrossValidator.__partition(nrow, self.folds)

        stats = PerfStats()
        lret = []
        xtra = []

        for f in xrange(self.folds):
            inpart = [i for i in xrange(nrow) if p[i] != f]
            outpart = [i for i in xrange(nrow) if p[i] == f]

            xin = x[inpart, :]
            yin = y[inpart]

            xout = x[outpart, :]
            yout = y[outpart]

            # print 'in:', xin.shape[0], 'out:', xout.shape[0], 'kwargs:', kwargs

            classifier = self.classifiercls(**kwargs)

            l = apply(getattr(classifier, self.__learnfunc), (xin, yin))
            if l is not None:
                lret.append(l)

            preds = apply(getattr(classifier, self.__predictfunc), (xout,))

            # do this after both learning and prediction just in case either performs some necessary computation
            if extra is not None:
                if isinstance(extra, str):
                    xtra.append(getattr(classifier, extra)())
                elif isinstance(extra, FunctionType):
                    xtra.append(extra(classifier))

            stats.append(*ystoconfusionmatrix(yout, preds))

        return {
            'learn': lret if len(lret) else None,
            'stats': stats,
            'extra': xtra if len(xtra) else None
        }
