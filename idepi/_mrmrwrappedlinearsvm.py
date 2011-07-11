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

from _mrmr import DiscreteMrmr
from _linearsvm import LinearSvm


__all__ = ['MrmrWrappedLinearSvm']


class MrmrWrappedLinearSvm(object):

    def __init__(self, num_features=10, method=DiscreteMrmr.MID, svm_type='c_svc',
                 degree=3, gamma=0.001, coef0=0, C=1, nu=0.5, eps=0.001, p=0.1,
                 cache_size=100, shrinking=True, probability=False, weight={}):
        self.__mrmr = Mrmr(num_features, method)
        self.__lsvm = LinearSvm(svm_type, degree, gamma, coef0, C, nu, eps, p,
                                cache_size, shrinking, probability, weight)

    def learn(self, x, y):
        self.__mrmr.select(x, y)
        x = self.__mrmr.subset(x)
        self.__lsvm.learn(x, y)

    def predict(self, x):
        x = self.__mrmr.subset(x)
        return self.__lsvm.predict(x)

    def subset(self, x):
        # raise ValueError('You need to learn() a model before you can reduce() other datasets.')
        return self.__mrmr.subset(x)

    def features(self):
        # raise ValueError('You need to learn() a model before you can determine the selected features().')
        return self.__mrmr.features()

    def weights(self):
        return self.__lsvm.weights()
