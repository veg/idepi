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

from warnings import warn

import numpy as np


__all__ = ['WrappedRegressor']


class WrappedRegressor(object):

    def __init__(self, selectorcls, regressorcls, **kwargs):

        kwargs0 = {}
        kwargs1 = {}
        for k, v in kwargs.items():
            if k[-1] == '0':
                kwargs0[k[:-1]] = v
            else:
                if k[-1] == '1':
                    k = k[:-1]
                kwargs1[k] = v

        self.__selector = selectorcls(**kwargs0)
        self.__regressor = regressorcls(**kwargs1)
        self.__beta0 = None
        self.__beta = None
        self.__idxs = None

    def learn(self, x, y):
        if not isinstance(x, np.ndarray):
            raise ValueError("x must be an numpy 2d array")

        if not isinstance(y, np.ndarray):
            raise ValueError("y must be an numpy 2d array")

        if x.ndim > 2:
            raise ValueError("x must be an 2d array")

        if x.shape[0] != y.shape[0]:
            raise ValueError("x and y are not aligned")

        self.__selector.learn(x, y)
        self.__beta0 = self.__selector.beta()

        self.__idxs = [0]
        for i in xrange(1, len(self.__beta0)):
            if self.__beta0[i] != 0.0:
                self.__idxs.append(i)
        x = x[:, self.__idxs]

        self.__regressor.learn(x, y)

        # reconstruct the full beta
#         beta = np.zeros(len(self.__beta0))
#         beta_ = self.__regressor.beta()
#         for i in xrange(len(self.__idxs)):
#             beta[self.__idxs[i]] = beta_[i]
        self.__beta = self.__regressor.beta()

    def beta0(self):
        return np.copy(self.__beta0)

    def beta(self):
        return np.copy(self.__beta)

    def pred(self, x):
        x = x[:, self.__idxs]
        return self.__regressor.pred(x)

    def selected(self):
        return np.copy(self.__idxs)

    def steps(self):
        if 'steps' in dir(self.__regressor):
            return self.__regressor.steps()
        warn('%s has no method `steps()\'' % self.__regressor.__class__.__name__)
        return None
