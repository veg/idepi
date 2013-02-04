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

from collections import Iterable
from math import sqrt


__all__ = ['NormalValue']


class NormalValue(list):

    def __init__(self, dtype, L=[], name=None):
        assert(isinstance(L, Iterable))
        assert(all(isinstance(l, dtype) for l in L))
        super(NormalValue, self).__init__(L)
        self.__name = name
        self.__dtype = dtype
        self.__mean = 0.0
        self.__std = 0.0
        self.__valid = False

    @property
    def name(self):
        return self.__name

    @property
    def mean(self):
        NormalValue.__compute(self)
        return self.__mean

    @property
    def std(self):
        NormalValue.__compute(self)
        return self.__std

    def __iadd__(self, values):
        NormalValue.extend(self, values)
        return self

    def __imul__(self, value):
        assert(isinstance(value, self.__dtype))
        self.__valid = False
        for i in range(len(self)):
            self[i] *= value
        return self

    def __mul__(self, value):
        assert(isinstance(value, self.__dtype))
        return NormalValue(self.__dtype, (v * value for v in self), self.name if hasattr(self, 'name') else None)

    def add(self, value):
        NormalValue.append(self, value)

    def append(self, value):
        assert(isinstance(value, self.__dtype))
        self.__valid = False
        super(NormalValue, self).append(value)

    def extend(self, values):
        assert(isinstance(values, list))
        assert(all(isinstance(v, self.__dtype) for v in values))
        self.__valid = False
        super(NormalValue, self).extend(values)

    def __compute(self):
        # keep the comprehensions or infinite recursions
        mean_ = (
            sum(self) / len(self)
            if len(self) else 0
            )
        std_ = 0.0
        for v in self:
            std_ += (mean_ - v) ** 2
        std_ = sqrt(std_)
        self.__mean = mean_
        self.__std = std_
        self.__valid = True

    def __ge__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mean >= other.mean else False

    def __gt__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mean > other.mean else False

    def __le__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mean <= other.mean else False

    def __lt__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mean < other.mean else False

    def __eq__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mean == other.mean and self.std == other.std else False

    def __repr__(self):
        return str((self.mean, self.std))

    def __str__(self):
        return '%g +- %g' % (self.mean, self.std)

    def __unicode__(self):
        return NormalValue.sprintf(self)

    def sprintf(self, fmt='%g \xb1 %g'):
        return fmt % (self.mean, self.std)
