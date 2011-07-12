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

from exceptions import ValueError
from itertools import chain
from math import ceil, floor
from random import randint, shuffle
from sys import exit, stderr
from types import ListType, IntType

import numpy as np

from Bio import AlignIO


__all__ = ['SeqTable']


class SeqTable(object):

    def __init__(self, alignment, alphabet, keep_func=None, skip_func=None):

        self.__alph = alphabet
        self.__alignment = alignment

        self.__keep = SeqTable.__filter(self, keep_func)
        self.__skip = SeqTable.__filter(self, skip_func)
        self.__rest = sorted(set(xrange(len(self.__alignment))) - (self.__keep | self.__skip))
        
        self.nrow = len(self.__alignment) - len(self.__keep) - len(self.__skip)
        self.ncol = self.__alignment.get_alignment_length()
        
        self.__fulldata = np.zeros((self.nrow, self.ncol), dtype=self.__alph.todtype)
        self.data = self.__fulldata
        self.cols = np.zeros((self.ncol, len(self.__alph)), dtype=int)
        self.cols[:, :] = np.sum(self.data, axis=0)

        SeqTable.__fill_data(self.__alignment, self.__alph, self.__fulldata, self.__keep, self.__skip)

        self.__folds = 1
        self.__rowlabels = -np.ones((len(alignment),), dtype=int)

        self.__mask = None

    @staticmethod
    def __fill_data(alignment, alphabet, data, keep_idxs, skip_idxs):
        _, ncol, _ = data.shape
        for i in xrange(len(alignment)):
            if i in keep_idxs or i in skip_idxs:
                continue
            row = alignment[i].seq
            for j in xrange(ncol):
                v = row.seq[j]
                data[i, j, alphabet[v]] = True

    def __filter(self, func):
        idxs = set()
        if func is not None:
            if type(func) is ListType:
                for f in func:
                    for i in xrange(len(self.__alignment)):
                        if apply(f, (self.__alignment[i].id,)):
                            idxs.add(i)
            else:
                for i in xrange(len(self.__alignment)):
                    if apply(func, (self.__alignment[i].id,)):
                        idxs.add(i)
        return idxs

    @staticmethod
    def __partition(l, folds):
        npf = int(floor(l / folds)) # num per fold
        r = l % folds
        p = list(chain(*[[i] * npf for i in xrange(folds)])) + range(r)
        shuffle(p)
        assert(len(p) == l)
        return p

    def partition(self, folds):
        if folds > self.nrow:
            raise ValueError('Unable to create more partitions than available rows')
        self.__rowlabels[self.__rest] = SeqTable.__partition(self.nrow, folds)

    def set_row_fold(self, i, j):
        self.__rowlabels[i] = j

    @property
    def alignment(self):
        if self.__mask is None or len(self.__mask) == 0:
            # this is tricky, so try to follow: keep everything not in skip indices, but add back in things in keep_indices
            return [self.__alignment[i] for i in self.__rest] + [self.__alignment[i] for i in sorted(self.__keep)]
            # this one is even trickier: keep everything that isn't in the current mask or skip indices, but add back in the things in the keep indices
        return [self.__alignment[i] for i in self.__rest if self.__rowlabels[i] not in self.__mask] + \
               [self.__alignment[i] for i in sorted(self.__keep)]

    def mask(self, i):
        if type(i) is int:
            self.__mask = set([i])
        elif type(i) is list:
            self.__mask = set(i)
        elif type(i) is set:
            self.__mask = i
        elif i is None:
            self.__mask = None
        else:
            raise ValueError('Mask must be an int, list of ints, set of ints, or None')
        if self.__mask is None or len(self.__mask) == 0:
            self.data = self.__fulldata
        else:
            self.data = self.__fulldata[sorted(self.__mask)]
        self.cols[:, :] = np.sum(self.data, axis=0)

    def unmask(self):
        self.__mask = None
        self.data = self.__fulldata
        self.cols[:, :] = np.sum(self.data, axis=0)
