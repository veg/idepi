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

from exceptions import AssertionError
from itertools import chain

import numpy as np


__all__ = ['Regressor']


class Regressor(object):

    def __init__(self, data, regressorcls=RidgeLar, *args, **kwargs):

        self.__method = method
        self.__regressor = regressorcls(*args, **kwargs)
        self.__learned = False

        self.__feature_names = data.feature_names
        self.__x, self.__y = Regressor.__smldata_to_xy(data)
        self.__x_active, self.__y_active = self.__x.copy(), self.__y.copy()
        self.__x_inactive, self.__y_inactive = None, None

        # define the normalization constants
        # define these here to make pylint be quiet
        self.__xbar = None
        self.__xvar = None
        self.__ybar = None

        self.__folds = None
        self.__partition = None

    def __normalize(self):
        nperr = np.geterr()
        np.seterr(divide='ignore')
        
        try:
            if self.__method not in ('gd',):
                __zero = pow(10.0, -9.0) # numerically close enough to 0.0
                ncol = self.__x_active.shape[1]
    
                # normalize y and validate
                self.__ybar = np.mean(self.__y_active)
                self.__y_active -= self.__ybar
                assert(np.abs(np.mean(self.__y_active)) < __zero)
    
                # normalize x and validate
                self.__xbar = np.zeros(ncol, dtype=float)
                self.__xvar = np.zeros(ncol, dtype=float)
                for j in xrange(ncol):
                    self.__xbar[j] = np.mean(self.__x_active[:, j])
                    self.__x_active[:, j] -= self.__xbar[j]
                    assert(np.abs(np.mean(self.__x_active[:, j])) < __zero)
                    if self.__method not in ('ridge', 'gd'):
                        self.__xvar[j] = np.sqrt(sum(pow(self.__x_active[:, j], 2.0)))
                        if self.__xvar[j] != 0.0:
                            self.__x_active[:, j] /= self.__xvar[j]
                            try:
                                assert(np.abs(sum([pow(i, 2.0) for i in self.__x_active[:, j]]) - 1.0) < __zero)
                            except AssertionError, e:
                                print u'\u03c3: %.4g, \u03a3x\u00b2: %.4g' % (self.__xvar[j], sum([pow(i, 2.0) for i in self.__x_active[:, j]]))
        finally:
            np.seterr(**nperr)

    def partition(self, l, folds):
        npf = int(np.floor(l / folds)) # num per fold
        r = len % folds
        p = list(chain(*[[i] * npf for i in xrange(folds)])) + range(r)
        np.random.shuffle(p)
        assert(len(p) == self.__x.shape[0])
        self.__folds = folds
        self.__partition = p

    def mask(self, fold):
        assert(0 <= fold and fold < self.__folds)
        active = [i for i in xrange(self.__x.shape[0]) if self.__partition[i] != fold]  
        inactive = [i for i in xrange(self.__x.shape[0]) if self.__partition[i] == fold]
        # don't need to issue a copy() here
        self.__x_active = self.__x[active, :]
        self.__y_active = self.__y[active]
        self.__x_inactive = self.__x[inactive, :]
        self.__y_inactive = self.__y[inactive]

    def unmask(self):
        self.__x_active, self.__y_active = self.__x.copy(), self.__y.copy()
        self.__x_inactive, self.__y_inactive = None, None

    @staticmethod
    def __smldata_to_xy(data):
        ncol = len(data.feature_names)

        x = []
        y = []

        for row in data:
            x.append([row.features[j] if j in row.features else 0. for j in xrange(ncol)])
            y.append(row.value)

        x = np.array(x).astype(float)
        x = np.hstack((np.ones((x.shape[0], 1), dtype=float), x)) # add a constant term (y-intercept) as the 0th column
        y = np.array(y).astype(float)

        return x, y

    def __normalize_for_predict(self, x, y):
        if self.__xbar is None or self.__xvar is None or self.__ybar is None:
            raise ValueError('normalization constants (xbar, xvar, yvar) are unset')

        nperr = np.geterr()
        np.seterr(divide='ignore')

        ncol = x.shape[1]

        if self.__method not in ('gd',):
            y -= self.__ybar # LAR, LASSO, ElasticNet: responses have mean 0
            # skip index 0 because it's our y-intercept
            for j in xrange(1, ncol):
                x[:, j] -= self.__xbar[j] # LAR, LASSO, ElasticNet: covariates have mean 0
                if self.__method not in ('ridge', 'gd'):
                    if self.__xvar[j] != 0.0:
                        x[:, j] /= self.__xvar[j] # LAR, LASSO: covariates have unit length

        np.seterr(**nperr)

        return x, y

    @property
    def intercept(self):
        if not self.__learned:
            Regressor.learn(self)
        return self.__regressor.beta()[0]

    @property
    def weights(self):
        if not self.__learned:
            Regressor.learn(self)
        beta = self.__regressor.beta()
        # ignore the y-intercept
        return dict(zip(self.__feature_names, beta[1:]))

    @property
    def selected(self):
        if not self.__learned:
            Regressor.learn(self)
        selected = self.__regressor.selected()
        # -1 refers to the y-intercept
        return np.array([s - 1 for s in selected])

    def learn(self):
        self.__regressor.learn(self.__x_active, self.__y_active)
        self.__learned = True

    def predict(self, data):
        x, y = Regressor.__smldata_to_xy(data)
        x, y = Regressor.__normalize_for_predict(self, x, y)
        if not self.__learned:
            Regressor.learn(self)
        return self.__regressor.pred(x)

    def test(self, data):
        nperr = np.geterr()
        np.seterr(divide='ignore')

        if not self.__learned:
            Regressor.learn(self)
        x, y = Regressor.__smldata_to_xy(data)
        x, y = Regressor.__normalize_for_predict(self, x, y)
        yhat = self.__regressor.pred(x)
        sse = sum(pow(y - yhat, 2.0))
        ybar = np.mean(y)
        sst = sum(pow(y - ybar, 2.0))
        r2 = 1.0 - (sse / sst)
        nless1 = len(y) - 1
        p = len([1 for i in self.weights.values() if i != 0.0]) - 1 # - 1 to avoid counting the constant term
        mse = sse / (nless1 - p) # count the the full N
        rmse = np.sqrt(mse)
        rbar2 = 1.0 - (1.0 - r2) * nless1 / (nless1 - p)

        np.seterr(**nperr)

        return {
            u'R\u0304\u00b2   ': rbar2,
            u'R\u00b2   ': r2,
            # u'SSE  ': sse,
            # u'MSE  ': mse,
            u'RMSE ': rmse,
        }
