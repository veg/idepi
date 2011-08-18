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
from os import close, remove
from os.path import exists
from re import sub
from shutil import copyfile
from subprocess import Popen, PIPE
from sys import stderr
from tempfile import mkstemp

import numpy as np

from mlpy import LibSvm


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
        fd, self.__modelfile = mkstemp(); close(fd)
        self.__lsvm.save_model(self.__modelfile)
        model = LinearSvmModel(self.__modelfile)
        remove(self.__modelfile)
        return model.weights()

    def __del__(self):
        if self.__modelfile is not None:
            if exists(self.__modelfile):
                remove(self.__modelfile)


class LinearSvmModel(object):

    def __init__(self, modelfile, feature_names=None):
        self.__feature_weights = {}
        self.__modelfile = modelfile
        self.__feature_names = feature_names
        self.type = None
        fh = open(self.__modelfile, 'rU')
        mode = None
        f = 0 # for LIBLINEAR models
        for line in fh:
            line = line.strip().upper()

            if line[:11] == 'KERNEL_TYPE':
                assert(line.split(' ', 1)[1] == 'LINEAR')

            # get the mode from the line right before the SVs (LIBSVM) or weight-vector (LIBLINEAR)
            if line == 'SV' or line == 'W':
                mode = line
                if mode == 'SV':
                    self.type = 'LIBSVM'
                elif mode == 'W':
                    self.type = 'LIBLINEAR'
                continue

            # if we don't have a mode yet or the line is empty skip
            if mode is None or line == '':
                continue

            vals = line.split()
            try:
                weight = float(vals[0])
            except ValueError:
                print >> stderr, 'ERROR: broken SVM model, skipping line \'%s\'' % line
                with open(self.__modelfile) as fh:
                    print >> stderr, fh.read()
                raise ValueError

            if mode == 'SV':
                if len(vals) <= 1:
                    # empty vector is useless, so keep going
                    continue
                for v in vals[1:]:
                    f, w = v.split(':')
                    # subtract 1 from here because model indices are 1-indexed for LIBSVM format
                    f = int(f)-1
                    if not f in self.__feature_weights:
                        self.__feature_weights[f] = 0.
                    self.__feature_weights[f] += weight * float(w)

            elif mode == 'W':
                if len(vals) > 1:
                    weight += sum([float(v) for v in vals[1:]])
                self.__feature_weights[f] = weight
                f += 1

    def weights(self):
        features = np.zeros((max(self.__feature_weights.keys()) + 1,), dtype=float)
        for k, v in self.__feature_weights.items():
            features[k] += v
        return features
