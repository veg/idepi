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

from operator import itemgetter
from os import remove
from os.path import exists
from re import sub
from shutil import copyfile
from subprocess import Popen, PIPE
from sys import stderr
from tempfile import mkstemp

from mlpy import LibSvm
from numpy import mean, std

from _smldata import SmlData
from _linearsvmmodel import LinearSvmModel


__all__ = ['LinearSvm']


class LinearSvm(object):

    def __init__(self, svm_type='c_svc', degree=3, gamma=0.001, coef0=0,
                 C=1, nu=0.5, eps=0.001, p=0.1, cache_size=100, shrinking=True,
                 probability=False, weight={}):
        self.__lsvm = LibSvm(svm_type, 'linear', degree, gamma, coef0, C, nu,
                             eps, p, cache_size, shrinking, probability,
                             weight)
        self.__modelfile = None
        self.__computed = False

    def learn(self, x, y):
        assert(x.shape[1] == 10)
        ret = self.__lsvm.learn(x.astype(float), y.astype(float))
        self.__computed = True
        if ret is not None:
            return ret
   
    def predict(self, x):
        return self.__lsvm.pred(x.astype(float))

    def predict_probabilities(self, x):
        return self.__lsvm.pred_probability(x.astype(float))

    def weights(self):
        if self.__computed == False:
            raise StandardError('No SVM model computed')
        self.__modelfile = mkstemp()[1]
        self.__lsvm.save_model(self.__modelfile)
        model = LinearSvmModel(self.__modelfile)
        remove(self.__modelfile)
        return model.weights()

    def __del__(self):
        if self.__modelfile is not None:
            if exists(self.__modelfile):
                remove(self.__modelfile)
