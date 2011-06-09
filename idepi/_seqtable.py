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
from math import ceil
from random import randint
from sys import exit, stderr
from Bio import AlignIO
from types import ListType, IntType


__all__ = ['SeqTable']


class Column(object):

    def __init__(self, n):
        self._values = dict()
        self.number = n

    def __contains__(self, key):
        if key in self._values:
            return True
        else:
            return False

    def __delitem__(self, key):
        del self._values[key]

    def __getitem__(self, key):
        return self._values[key]

    def __len__(self):
        return len(self._values.keys())

    def __setitem__(self, key, value):
        self._values[key] = value

    def values(self):
        return self._values.keys()

    def counts(self):
        return dict(self._values)

    def count(self):
        return sum(self._values.values())


class SeqTable(object):
    DNA_ALPHABET   = 'ACGTUN-'
    AMINO_ALPHABET = 'ACGILMPSTVDENQFWYHKRX-'

    # underscore must go first or re.compile blows up
    SPACE = '_.-='

    def __init__(self, filename, alphabet, keep_func = None, skip_func = None):
        if not alphabet in (self.DNA_ALPHABET, self.AMINO_ALPHABET):
            print >> stderr, 'ERROR: alphabet must be one of SeqTable.DNA_ALPHABET or SeqTable.AMINO_ALPHABET'
            raise ValueError

        self._fh = open(filename, 'rU')
        self.alphabet = alphabet
        self.alignment = AlignIO.read(self._fh, 'stockholm')

        self.keep_indices = SeqTable.__filter(self, keep_func)
        self.skip_indices = SeqTable.__filter(self, skip_func)

        self.num_columns = self.alignment.get_alignment_length()
        self.num_rows = len(self.alignment)

        self.num_partitions = 1
        self.row_partitions = [-1] * self.num_rows

        self.columns = {}
        self.curr_mask = None
        self.loaded = False

    def __filter(self, func):
        indices = []
        if func:
            if type(func) is ListType:
                for f in func:
                    for i in xrange(len(self.alignment)):
                        if f(self.alignment[i].id):
                            indices.append(i)
            else:
                for i in xrange(len(self.alignment)):
                    if func(self.alignment[i].id):
                        indices.append(i)
        return indices

    def fill_columns(self):
        if self.loaded:
            return

        self.columns.clear()

        for r in self.alignment:
            for p in r.seq:
                if p in self.SPACE and p != '-':
                    print >> stderr, 'ERROR: space character not matching '-' found!'
                    raise ValueError

            for i in xrange(self.alignment.get_alignment_length()):
                v = r.seq[i]

                if v not in self.alphabet:
                    v = 'X'

                if not i in self.columns:
                    self.columns[i] = Column(i)

                # initialize if not already there
                if not v in self.columns[i]:
                    self.columns[i][v] = 0

                self.columns[i][v] += 1

        self.loaded = True

    def partition(self, num_partitions):
        self.num_partitions = num_partitions

        # set div to div++ to accomodate remainder logic below
        nr_rows = self.num_rows - len(self.skip_indices) - len(self.keep_indices)
        div = int(ceil(nr_rows / self.num_partitions))
        rem = int(nr_rows % self.num_partitions)

        indices = [i for i in xrange(self.num_rows) if i not in (self.keep_indices + self.skip_indices)]

        if len(indices) > num_partitions:
            for i in xrange(num_partitions):
                for j in xrange(div):
                    try:
                        k = randint(0, len(indices)-1)
                    except ValueError:
                        k = 0
                    try:
                        self.row_partitions[indices[k]] = i
                        del indices[k]
                    except IndexError, e:
                        print >> stderr, 'ERROR: ran out of indices before we should have!'
                        raise e
                # keep decrementing remainder until we hit 0, then drop divisor by 1 and keep going
                rem -= 1
                if rem == 0:
                    div -= 1
        elif len(indices) == num_partitions:
            c = 0
            for i in indices:
                self.row_partitions[i] = c
                c += 1
        else:
            print >> stderr, 'ERROR: can\'t create more partitions than available rows'
            exit(-1)

    def set_row_fold(self, i, j):
        self.row_partitions[i] = j

    def rows(self):
        if not self.loaded:
            SeqTable.fill_columns(self)
        if self.curr_mask is None or self.curr_mask == []:
            # this is tricky, so try to follow: keep everything not in skip indices, but add back in things in keep_indices
            return [self.alignment[i] for i in xrange(len(self.alignment)) if i not in (self.keep_indices + self.skip_indices)] + [self.alignment[i] for i in self.keep_indices]
            # this one is even trickier: keep everything that isn't in the current mask or skip indices, but add back in the things in keep_indices
        return [self.alignment[i] for i in xrange(len(self.alignment)) if self.row_partitions[i] not in self.curr_mask and i not in (self.keep_indices + self.skip_indices)] + \
               [self.alignment[i] for i in self.keep_indices]

    def mask(self, i):
        if type(i) is IntType:
            self.curr_mask = [i]
        elif type(i) is ListType:
            self.curr_mask = i
        else:
            print >> stderr, 'ERROR: mask must be of type int or a list of ints'
            raise ValueError
        self.loaded = False

    def unmask(self):
        self.curr_mask = None
        self.loaded = False

    def __del__(self):
        self._fh.close()
