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

from sys import stderr

import numpy as np


__all__ = ['SmlData']


class Record:

    def __init__(self, value, features):
        self.value = value
        self.features = dict(features)


class SmlData:

    def __init__(self, feature_names):
        self.__data = []
        self.feature_names = feature_names

    def add(self, value, features):
        if isinstance(value, list) and isinstance(features, list):
            if len(value) == len(features) and len(features) != 0:
                self.__data.extend([Record(value[i], features[i]) for i in range(len(value))])
            elif len(value) != len(features):
                print('ERROR: number of values and features don\'t match for multiple-add', file=stderr)
                exit(-1)
        elif not isinstance(features, dict):
            print('ERROR: feature argument takes a dict()', file=stderr)
            exit(-1)
        else:
            self.__data.append(Record(value, features))

    def __contains__(self, key):
        if key >= 0 and len(self.__data) > key:
            return True
        else:
            return False

    def __delitem__(self, key):
        if SmlData.__contains__(self, key):
            del self.__data[key]
        else:
            raise IndexError

    def __getitem__(self, key):
        if SmlData.__contains__(self, key):
            return self.__data[key]
        else:
            raise IndexError

    def __iter__(self):
        return self.__data.__iter__()

    def __len__(self):
        return len(self.__data)

    def __setitem__(self, key, value):
        if SmlData.__contains__(self, key):
            self.__data[key] = value
        else:
            raise IndexError

    def save_tab(self, filename, target):
        fh = open(filename, 'w')
        print('%s\t%s' % ('\t'.join([self.feature_names[i] for i in range(len(self.feature_names))]), target), file=fh)
        print('%s\td' % '\t'.join(['d' for i in range(len(self.feature_names))]), file=fh)
        print('%s\tclass' % '\t'.join(['' for i in range(len(self.feature_names))]), file=fh)
        for r in self:
            line = '\t'.join(['1' if i in r.features else '0' for i in range(len(self.feature_names))])
            if r.value is None:
                raise ValueError('ERROR: save_tab doesn\'t yet support value-less SmlData')
            print('%s\t%d' % (line, r.value), file=fh)
        fh.close()

    def tondarrays(self):
        y = np.array([r.value for r in self])
        x = np.array([[1 if i in r.features else 0 for i in range(len(self.feature_names))] for r in self])

        return x, y
