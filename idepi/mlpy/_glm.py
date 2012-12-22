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

import numpy as np


__all__ = ['Glm']


class Glm:

    __linkfuns = {
        'id': lambda x: x,
        'inv': lambda x: 1.0 / x,
        'invsq': lambda x: 1.0 / np.square(x),
        'log': lambda x: np.log(x),
        'logit': lambda x: np.log(x / (1.0 - x)),
    }

    __invlinkfuns = {
        'id': lambda x: x,
        'inv': lambda x: 1.0 / x,
        'invsq': lambda x: 1.0 / np.sqrt(x),
        'log': lambda x: np.exp(x),
        'logit': lambda x: 1.0 / (1.0 + np.exp(-x)),
    }

    __families = {
        'gaussian': 'id',
        'normal': 'id',
        'exponential': 'inv',
        'gamma': 'inv',
        'invgauss': 'invsq',
        'poisson': 'log',
        'binomial': 'logit',
        'multinomial': 'logit',
    }

    # TODO: finish these initializer functions
    __initfuns = {
        'gaussian': lambda x: x,
        'normal': lambda x: x,
        'exponential': lambda x: x, # this is probably x, but I don't know
        'gamma': lambda x: x,
        'invgauss': lambda x: x,
        'poisson': lambda x: x + 0.1,
        'binomial': lambda x: (x + 0.5) / 2.0,
        'multinomial': lambda x: (x + 0.5) / 2.0, # this is probably correct, but again I don't know
    }

    __varfuns = {
        'gaussian': lambda x: np.ones(x.shape, dtype=float),
        'normal': lambda x: np.ones(x.shape, dtype=float),
        'exponential': lambda x: 1.0 / pow(x, 2),
        'gamma': lambda x: pow(x, 2),
        'invgauss': lambda x: pow(x, 3),
        'poisson': lambda x: x,
        'binomial': lambda x: x * (1.0 - x),
        'multinomial': lambda x: x * (1.0 - x),
    }

    def __init__(self, dist='normal'):

        if dist not in self.__families.keys():
            raise ValueError('dist `%s\' must be one of (%s)' % (dist, ', '.join(self.__families.keys())))

        self.dist = dist
        self.name = self.__families[dist]

        self.link = self.__linkfuns[self.name]
        self.invlink = self.__invlinkfuns[self.name]
        self.initialize = self.__initfuns[self.dist]
        self.var = self.__varfuns[self.dist]
