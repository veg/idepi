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

from mlpy import Lar

from _dantzig import Dantzig
from _wrappedregressor import WrappedRegressor


__all__ = ['LarDantzig']


class LarDantzig(WrappedRegressor):
    def __init__(self, tol=0.001, lam=1.0, m=None, **kwargs):
        if m is None:
            raise ValueError('you must specify a maximum number of iterations m')
        if 'tol0' in kwargs:
            lam=kwargs['tol0']
            del kwargs['tol0']
        if 'lam0' in kwargs:
            lam=kwargs['lam0']
            del kwargs['lam0']
        super(LarDantzig, self).__init__(selectorcls=Dantzig, regressorcls=Lar, tol0=tol, lam0=lam, m1=m, **kwargs)
