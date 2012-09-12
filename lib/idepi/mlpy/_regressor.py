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

from __future__ import division, print_function

import numpy as np

from mlpy import ElasticNet, LARS, Ridge

from ._dantzig import Dantzig
from ._doubledantzig import DoubleDantzig
from ._lardantzig import LarDantzig
from ._linearsvr import LinearSvr
from ._ridgedantzig import RidgeDantzig
from ._ridgelar import RidgeLar


__all__ = ['Regressor', 'regressor_classes']


class Regressor(object):
    DANTZIG = Dantzig
    DBLDANTZIG = DoubleDantzig
    ELASTICNET = ElasticNet
    LAR = LARS
    LARDANTZIG = LarDantzig
    LINEARSVR = LinearSvr
    RIDGE = Ridge
    RIDGEDANTZIG = RidgeDantzig
    RIDGELAR = RidgeLar

    def __init__(self, regressorcls=RidgeLar, logspace=False, *args, **kwargs):

        self.__regressor = regressorcls(*args, **kwargs)
        self.__learned = False
        self.__logspace = logspace

        assert(np.log(np.exp(1)) == 1.)

        # define the normalization constants
        # define these here to make pylint be quiet
        self.__xbar = None
        self.__xvar = None
        self.__ybar = None

    @staticmethod
    def __normalize(x, y):
        ZERO = pow(10.0, -9.0) # numerically close enough to 0.0

        nperr = np.geterr()
        np.seterr(divide='ignore')

        x = Regressor.__addintercept(x)

        # shouldn't be necessary
        # y = y.astype(float)

        try:
#           if self.__method not in ('gd',):
            ncol = x.shape[1]

            # normalize y and validate
            ybar = np.mean(y)
            y -= ybar
            assert(np.abs(np.mean(y)) < ZERO)

            # normalize x and validate
            xbar = np.zeros(ncol, dtype=float)
            xvar = np.zeros(ncol, dtype=float)
            for j in range(1, ncol):
                xbar[j] = np.mean(x[:, j])
                x[:, j] -= xbar[j]
                assert(np.abs(np.mean(x[:, j])) < ZERO)
#               if self.__method not in ('ridge', 'gd'):
                xvar[j] = np.sqrt(sum(pow(x[:, j], 2.0)))
                if xvar[j] != 0.0:
                    x[:, j] /= xvar[j]
                    try:
                        assert(np.abs(sum([pow(i, 2.0) for i in x[:, j]]) - 1.0) < ZERO)
                    except AssertionError:
                        print('\u03c3: %.4g, \u03a3x\u00b2: %.4g' % (xvar[j], sum([pow(i, 2.0) for i in x[:, j]])))
        finally:
            np.seterr(**nperr)

        return x, y, xbar, xvar, ybar

    @staticmethod
    def __addintercept(x):
        nrow, ncol = x.shape
        ncol += 1
        newx = np.empty((nrow, ncol), dtype=float)
        newx[:, 0 ] = 1
        newx[:, 1:] = x
        return newx

    @staticmethod
    def __normalize_for_predict(x, xbar, xvar):
        nperr = np.seterr(divide='ignore')

        x = Regressor.__addintercept(x)

        ncol = x.shape[1]

#        if self.__method not in ('gd',):
#        y -= self.__ybar # LARS, ElasticNet: responses have mean 0
        # skip index 0 because it's our y-intercept
        for j in range(1, ncol):
            x[:, j] -= xbar[j] # LARS, ElasticNet: covariates have mean 0
#            if self.__method not in ('ridge', 'gd'):
            if xvar[j] != 0.0:
                x[:, j] /= xvar[j] # LARS: covariates have unit length

        np.seterr(**nperr)

        return x

    def intercept(self):
        if not self.__learned:
            raise RuntimeError('No regression model computed')
        return self.__regressor.beta()[0]

    def weights(self):
        if not self.__learned:
            raise RuntimeError('No regression model computed')
        beta = self.__regressor.beta()
        # ignore the y-intercept
        return beta[1:]

    def selected(self):
        if not self.__learned:
            raise RuntimeError('No regression model computed')
        # -1 refers removes the y-intercept
        selected = self.__regressor.selected()
        return np.array([s - 1 for s in selected])

    def learn(self, x, y):
#         if np.any(y == 1.):
#             print >> stderr, 'Cannot perform logspace conversion, divide by zeros will result. Falling back to non-logspace regression.'
#             self.__logspace = False
        if self.__logspace:
            y = np.log(y)
        x, y, self.__xbar, self.__xvar, self.__ybar = Regressor.__normalize(x, y)
        self.__regressor.learn(x, y)
        self.__learned = True

    def predict(self, x):
        if not self.__learned:
            raise RuntimeError('No regression model computed')
        x = Regressor.__normalize_for_predict(x, self.__xbar, self.__xvar)
        y = self.__regressor.pred(x) + self.__ybar
        if self.__logspace:
            y = np.exp(y)
        return y

#     def test(self, data):
#         nperr = np.seterr(divide='ignore')
#
#         if not self.__learned:
#             Regressor.learn(self)
#         x, y = Regressor.__smldata_to_xy(data)
#         x, y = Regressor.__normalize_for_predict(self, x, y)
#         yhat = self.__regressor.pred(x)
#         sse = sum(pow(y - yhat, 2.0))
#         ybar = np.mean(y)
#         sst = sum(pow(y - ybar, 2.0))
#         r2 = 1.0 - (sse / sst)
#         nless1 = len(y) - 1
#         p = len([1 for i in self.weights.values() if i != 0.0]) - 1 # - 1 to avoid counting the constant term
#         mse = sse / (nless1 - p) # count the the full N
#         rmse = np.sqrt(mse)
#         rbar2 = 1.0 - (1.0 - r2) * nless1 / (nless1 - p)
#
#         np.seterr(**nperr)
#
#         return {
#             u'R\u0304\u00b2   ': rbar2,
#             u'R\u00b2   ': r2,
#             # u'SSE  ': sse,
#             # u'MSE  ': mse,
#             u'RMSE ': rmse,
#         }

regressor_classes = {
    'dantzig': Regressor.DANTZIG,
    'doubledantzig': Regressor.DBLDANTZIG,
    'elasticnet': Regressor.ELASTICNET,
    'lar': Regressor.LAR,
    'lardantzig': Regressor.LARDANTZIG,
    'linearsvr': Regressor.LINEARSVR,
    'ridge': Regressor.RIDGE,
    'ridgedantzig': Regressor.RIDGEDANTZIG,
    'ridgelar': Regressor.RIDGELAR,
}
