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
from numpy import mean, nan_to_num, var


__all__ = ['NormalValue']


class NormalValue(list):

    def __init__(self, dtype, L=[], name=None):
        assert(isinstance(L, Iterable))
        assert(all(isinstance(l, dtype) for l in L))
        super(NormalValue, self).__init__(L)
        self.__dtype = dtype
        if name is not None:
            self.name = name
        self.mu = 0.
        self.sigma = 0.
        NormalValue.__compute(self)

    def __iadd__(self, values):
        NormalValue.extend(self, values)
        return self

    def __imul__(self, value):
        assert(isinstance(value, self.__dtype))
        for i in range(len(self)):
            self[i] *= value
        NormalValue.__compute(self)
        return self

    def __mul__(self, value):
        assert(isinstance(value, self.__dtype))
        return NormalValue(self.__dtype, (v * value for v in self), self.name if hasattr(self, 'name') else None)

    def append(self, value):
        assert(isinstance(value, self.__dtype))
        super(NormalValue, self).append(value)
        NormalValue.__compute(self)

    def extend(self, values):
        assert(isinstance(values, list))
        assert(all(isinstance(v, self.__dtype) for v in values))
        super(NormalValue, self).extend(values)
        NormalValue.__compute(self)

    def __compute(self):
        if len(self) == 0:
            self.mu = 0.
            self.sigma = 0.
        else:
            self.mu = float(nan_to_num(mean(self)))
            self.sigma = float(nan_to_num(var(self)))

    def __ge__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mu >= other.mu else False

    def __gt__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mu > other.mu else False

    def __le__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mu <= other.mu else False

    def __lt__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mu < other.mu else False

    def __eq__(self, other):
        assert(isinstance(other, NormalValue))
        return True if self.mu == other.mu and self.sigma == other.sigma else False

    def __repr__(self):
        return str((self.mu, self.sigma))

    def __str__(self):
        return '%g +- %g' % (self.mu, sqrt(self.sigma))

    def __unicode__(self):
        return NormalValue.sprintf(self)

    def sprintf(self, format='%g \xb1 %g'):
        return format % (self.mu, sqrt(self.sigma))
