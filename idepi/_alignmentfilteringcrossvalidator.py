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

from _crossvalidator import CrossValidator
from _perfstats import PerfStats
from _util import ystoconfusionmatrix


__all__ = ['AlignmentFilteringCrossValidator']


# implement cross-validation interface here, grid-search optional
class AlignmentFilteringCrossValidator(CrossValidator):

    def __init__(self, classifiercls, filtercls, folds, extract_func, cv={}, af={}, ref_id_func=None, skip_func=None, discretize_func=None):
        if ref_id_func is None:
            ref_id_func = lambda x: False

        if skip_func is None:
            skip_func = lambda x: False

        self.__extract = extract_func
        self.__discretize = discretize_func

        self.__skip = skip_func
        self.__ref_id = ref_id_func

        self.filtercls = filtercls
        self.af = af

        super(AlignmentFilteringCrossValidator, self).__init__(classifiercls, folds, cv)

    def crossvalidate(self, alignment, cv={}, extra=None, af={}):
        if extra is not None:
            if not isinstance(extra, str) and not isinstance(extra, FunctionType):
                raise ValueError('the `extra\' argument takes either a string or a function.')

        af_kwargs = deepcopy(self.af)
        af_kwargs.update(af)

        cv_kwargs = deepcopy(self.cv)
        cv_kwargs.update(cv)

        ignore_idxs = set()
        refseq = None
        for i in xrange(len(alignment)):
            row = alignment[i]
            r, s = apply(self.__ref_id, (row.id,)), apply(self.__skip, (row.id,))
            if r or s:
                ignore_idxs.add(i)
                if r and refseq is None: 
                    refseq = row
                elif r:
                    raise RuntimeError('Refseq found twice!!!')

        # filter out all the shit we don't want
        alignment = [alignment[i] for i in xrange(len(alignment)) if i not in ignore_idxs]

        nrow = len(alignment)

        p = CrossValidator.__partition(nrow, self.folds)

        yextractor = ClassExtractor(
            self.__extract,
            lambda row: apply(self.__ref_id, (row.id,)) or apply(self.__skip, (row.id,)),
            self.__discretize
        )
        y = yextractor.extract(alignment)

        assert(len(y) == nrow)

        names = []
        stats = PerfStats()
        lret = []
        xtra = []

        for f in xrange(self.folds):
            infold = [alignment[i] for i in xrange(nrow) if p[i] != f]
            outfold = [alignment[i] for i in xrange(nrow) if p[i] == f]

            inpart = [i for i in xrange(nrow) if p[i] != f]
            outpart = [i for i in xrange(nrow) if p[i] == f]

            colfilter = self.filtercls(**af_kwargs)

            if refseq is not None:
                colnames = colfilter.learn(infold + [refseq]) 
            else:
                colnames = colfilter.learn(infold)
            
            if colnames is not None:
                names.append(colnames)

            xin = colfilter.filter(infold).astype(float)
            yin = y[inpart].astype(float)

            xout = colfilter.filter(outfold).astype(float)
            yout = y[outpart].astype(float)

            # print 'in:', xin.shape[0], 'out:', xout.shape[0], 'kwargs:', kwargs

            classifier = self.classifiercls(**cv_kwargs)

            l = getattr(classifier, self.__learnfunc)(xin, yin)
            if l is not None:
                lret.append(l)

            preds = getattr(classifier, self.__predictfunc)(xout)

            # do this after both learning and prediction just in case either performs some necessary computation
            if extra is not None:
                if isinstance(extra, str):
                    xtra.append(getattr(classifier, extra)())
                elif isinstance(extra, FunctionType):
                    xtra.append(extra(classifier))

            stats.append(*ystoconfusionmatrix(yout, preds))

        return {
            'names': names if len(names) else None,
            'learn': lret if len(lret) else None,
            'stats': stats,
            'extra': xtra if len(xtra) else None
        }
