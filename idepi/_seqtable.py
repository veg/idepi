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

import numpy as np

from Bio import AlignIO


__all__ = ['SeqTable']


class SeqTable(object):

    def __init__(self, alignment, alphabet, ref_id_func=None, skip_func=None):

        self.__alph = alphabet
        self.__alignment = alignment

        self.__ref_id = SeqTable.__filter(self, ref_id_func)
        self.__skip = SeqTable.__filter(self, skip_func)
        self.__keep = sorted(set(xrange(len(self.__alignment))) - (self.__ref_id | self.__skip))

        self.nrow = len(self.__alignment) - len(self.__ref_id) - len(self.__skip)
        self.ncol = self.__alignment.get_alignment_length()

        assert(len(self.__keep) == self.nrow)

        self.__fulldata = np.zeros((self.nrow, self.ncol), dtype=self.__alph.todtype())
        self.data = self.__fulldata
        self.cols = np.zeros((self.ncol, len(self.__alph)), dtype=int)

        SeqTable.__fill_data(self.__alignment, self.__alph, self.__fulldata, self.__ref_id, self.__skip)

        self.cols[:, :] = np.sum(self.data, axis=0)

        self.__folds = 1
        self.__rowlabels = -np.ones((len(alignment),), dtype=int)

        self.__mask = None

    @staticmethod
    def __fill_data(alignment, alphabet, data, keep_idxs, skip_idxs):
        _, ncol, _ = data.shape
        off = 0
        for i, row in enumerate(alignment):
            if i in keep_idxs or i in skip_idxs:
                off += 1
                continue
            seq = str(row.seq)
            for j in xrange(ncol):
                v = seq[j]
                data[i - off, j, alphabet[v]] = True

    def __filter(self, func):
        idxs = set()
        if func is not None:
            if isinstance(func, list):
                for f in func:
                    for i, row in enumerate(self.__alignment):
                        if apply(f, (row.id,)):
                            idxs.add(i)
            else:
                for i, row in enumerate(self.__alignment):
                    if apply(func, (row.id,)):
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
        self.__rowlabels[self.__keep] = SeqTable.__partition(self.nrow, folds)

    def set_row_fold(self, i, j):
        self.__rowlabels[i] = j

    @property
    def alignment(self):
        if self.__mask is None or len(self.__mask) == 0:
            # this is tricky, so try to follow: keep everything not in skip indices, but add back in the refseqs
            return [self.__alignment[i] for i in self.__keep] + [self.__alignment[i] for i in sorted(self.__ref_id)]
            # this one is even trickier: keep everything that isn't in the current mask or skip indices, but add back in the refseqs
        return [self.__alignment[i] for i in self.__keep if self.__rowlabels[i] not in self.__mask] + \
               [self.__alignment[i] for i in sorted(self.__ref_id)]

    def mask(self, i):
        if isinstance(i, int):
            self.__mask = set([i])
        elif isinstance(i, list):
            self.__mask = set(i)
        elif isinstance(i, set):
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
