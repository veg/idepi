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
from types import FunctionType

from _perfstats import *
from _normalvalue import NormalValue
from _crossvalidator import CrossValidator


__all__ = [
    'ACCURACY',
    'PPV',
    'PRECISION',
    'NPV',
    'SENSITIVITY',
    'RECALL',
    'SPECIFICITY',
    'TNR',
    'FSCORE',
    'MINSTAT',
    'GridSearcher'
]


# implement cross-validation interface here, grid-search optional
class GridSearcher(CrossValidator):

    def __init__(self, classifiercls, folds, cv={}, optstat=MINSTAT, gs={}):
        super(GridSearcher, self).__init__(classifiercls, folds, cv)
        self.gs = gs
        self.optstat = optstat
        self.classifier = None
        self.__computed = False

    def gridsearch(self, x, y, cv={}, extra=None):
        if extra is not None:
            if not isinstance(extra, str) and not isinstance(extra, FunctionType):
                raise ValueError('the `extra\' argument takes either a string or a function.') 

        ret = { 'stats': PerfStats() }

        if len(self.gs) == 1:
            k0, params = self.gs.items()[0]
            bestparams = {}
            for p0 in params:
                cv[k0] = p0
                r = GridSearcher.crossvalidate(self, x, y, cv=cv)
                if r['stats'].get(self.optstat) > ret['stats'].get(self.optstat):
                    ret = r
                    bestparams = deepcopy(cv)
            kwargs = deepcopy(self.cv)
            for k, v in bestparams.items():
                kwargs[k] = v
            ret['kwargs'] = kwargs
        elif len(self.gs) == 2:
            gsp = self.gs.items()
            k0, params0 = gsp[0]
            k1, params1 = gsp[1]
            bestparams = {}
            for p0 in params0:
                for p1 in params1:
                    cv[k0] = p0, cv[k1] = p1
                    r = GridSearcher.crossvalidate(self, x, y, cv=cv)
                    if r['stats'].get(self.optstat) > ret['stats'].get(self.optstat):
                        ret = r
                        bestparams = deepcopy(cv) 
            kwargs = deepcopy(self.cv)
            for k, v in bestparams.items():
                kwargs[k] = v
            ret['kwargs'] = kwargs 
        else:
            raise ValueError('We only support up to a 2D grid search at this time')

        ret['extra'] = GridSearcher.crossvalidate(self, x, y, cv=ret['kwargs'], extra=extra)['extra'] if extra is not None else None

#         print ret['kwargs']
#         print '\n'.join([str(s) for s in ret['stats'].tolist()])

        return ret

    def learn(self, x, y):
        gsret = GridSearcher.gridsearch(self, x, y)

        self.classifier = self.classifiercls(**gsret['kwargs'])
        # I don't like unmangling the private name, but here it is..
        lret = getattr(self.classifier, self._CrossValidator__learnfunc)(x, y)

        self.__computed = True

        return { 'gridsearch': gsret, 'learn': lret }

    def predict(self, x):
        if self.__computed == False:
            raise RuntimeError('No model computed')

        return getattr(self.classifier, self._CrossValidator__predictfunc)(x)
