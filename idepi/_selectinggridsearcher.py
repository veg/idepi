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

from sys import stderr
from copy import deepcopy

from _perfstats import PerfStats
from _gridsearcher import GridSearcher


__all__ = ['SelectingGridSearcher']


# implement cross-validation interface here, grid-search optional
class SelectingGridSearcher(GridSearcher):

    ACCURACY            = PerfStats.ACCURACY
    PPV, PRECISION      = PerfStats.PPV, PerfStats.PRECISION
    NPV                 = PerfStats.NPV
    SENSITIVITY, RECALL = PerfStats.SENSITIVITY, PerfStats.RECALL
    SPECIFICITY, TNR    = PerfStats.SPECIFICITY, PerfStats.TNR
    FSCORE              = PerfStats.FSCORE
    MINSTAT             = PerfStats.MINSTAT

    def __init__(self, classifiercls, selectorcls, folds, cv={}, mode=None, optstat=PerfStats.MINSTAT, gs={}, fs={}):
        super(SelectingGridSearcher, self).__init__(classifiercls, folds, cv, mode, optstat, gs)
        self.__selected = False
        self.selector = selectorcls(**fs)

    def select(self, x, y):
        self.selector.select(x, y)
        self.__selected = True

    def subset(self, x):
        if self.__selected == False:
            raise StandardError('Selection hasn\'t yet been performed.')
        return self.selector.subset(x)

    def features(self):
        if self.__selected == False:
            raise StandardError('Selection hasn\'t yet been performed.')
        return self.selector.features()

    def gridsearch(self, x, y, cv={}, extra=None):
        SelectingGridSearcher.select(self, x, y)
        x = self.selector.subset(x)
        return super(SelectingGridSearcher, self).gridsearch(x, y, cv=cv, extra=extra)

    def learn(self, x, y):
        SelectingGridSearcher.select(self, x, y)
        x = self.selector.subset(x)
        ret = super(SelectingGridSearcher, self).learn(x, y)
        if ret is not None:
            return ret

    def predict(self, x):
        if self.__selected == False:
            raise StandardError('You need to call learn() or select() before you can predict with this classifier.')
        x = self.selector.subset(x)
        return super(SelectingGridSearcher, self).predict(x)
