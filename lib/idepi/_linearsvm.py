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

from operator import itemgetter
from os import close, remove
from os.path import exists
from re import sub
from shutil import copyfile
from subprocess import Popen, PIPE
from tempfile import mkstemp
from warnings import warn

import numpy as np

from mlpy import LibSvm


__all__ = ['LinearSvm']


def bias_func(x):
    return np.hstack((x, np.ones((x.shape[0], 1), dtype=float)))

def nobias_func(x):
    return x.astype(float)


class LinearSvm(object):

    def __init__(self, svm_type='c_svc', degree=3, gamma=0.001, coef0=0,
                 C=1, nu=0.5, eps=0.001, p=0.1, cache_size=100, shrinking=True,
                 probability=False, weight={}, bias=False):
        self.__lsvm = LibSvm(svm_type, 'linear', degree, gamma, coef0, C, nu,
                             eps, p, cache_size, shrinking, probability,
                             weight)
        self.__bias = bias_func if bias else nobias_func
        self.__modelfile = None
        self.__computed = False

    def learn(self, x, y):
        ret = self.__lsvm.learn(self.__bias(x), y.astype(float))
        self.__computed = True
        if ret is not None:
            return ret

    def predict(self, x):
        return self.__lsvm.pred(self.__bias(x))

    def predict_probabilities(self, x):
        return self.__lsvm.pred_probability(self.__bias(x))

    def weights(self):
        if self.__computed == False:
            raise Exception('No SVM model computed')
        fd, self.__modelfile = mkstemp(); close(fd)
        self.__lsvm.save_model(self.__modelfile.encode('utf-8'))
        model = LinearSvmModel(self.__modelfile)
        remove(self.__modelfile)
        return model.weights()

    def __del__(self):
        if self.__modelfile is not None:
            if exists(self.__modelfile):
                try:
                    remove(self.__modelfile)
                except:
                    warn("file '%s' removed between existence check and call to unlink()" % self.__modelfile)


class LinearSvmModel(object):

    def __init__(self, modelfile, feature_names=None):
        self.__feature_weights = {}
        self.__modelfile = modelfile
        self.__feature_names = feature_names
        self.type = None
        fh = open(self.__modelfile, 'rU')
        mode = None
        rev = 1
        f = 0 # for LIBLINEAR models
        for line in fh:
            line = line.strip().upper()

            # if the labels are reversed, that means that libsvm switched the labelling internally, fix this
            if line[:5] == 'LABEL':
                arr = [float(v) for v in line[6:].split()]
                # if label 0 1 then things are backwards, so flip them
                if len(arr) == 2 and arr[0] < arr[1]:
                    rev = -1

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
                weight = float(vals[0]) * rev
            except ValueError:
                lines = None
                with open(self.__modelfile) as fh:
                    lines = fh.read().strip()
                raise ValueError("broken SVM model:  \n%s" % lines.replace('\n', '  \n'))

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
                    weight += sum(float(v) for v in vals[1:])
                self.__feature_weights[f] = weight
                f += 1

    def weights(self):
        numfeats = max(self.__feature_weights.keys()) + 1 if len(self.__feature_weights) else 0
        features = np.zeros((numfeats,), dtype=float)
        for k, v in self.__feature_weights.items():
            features[k] += v
        return features
