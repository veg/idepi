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

from warnings import warn

from cvxopt import matrix as co_matrix
from cvxopt import solvers as co_solvers

import numpy as np

from ._glm import Glm


__all__ = ['Dantzig']


class Dantzig:
    '''
    This code is basically a pythonic copy of
    gdscode.R provided by Gareth James (UCSC) as supplementary material
    to his paper with Radchenko, `A Generalized Dantzig Selector with
    Shrinkage Tuning', (2009). Biometrika, 98, 323-337.
    '''

    __i = -1 # set to -1 by default to make shit happen

    def __init__(self, tol=0.001, lam=0.001, m=None, dist='normal'):
        self.__tol = tol
        self.__lam = lam
        self.__beta = None
        self.__selected = None
        self.__glm = Glm(dist)
        self.__m = m
        self.__L = 0

    @classmethod
    def __linearoptim(cls, x, y, sigma):
        # save default numpy errors
        nperr = np.geterr()

        np.seterr(divide='ignore')

        # silence the LP solvers
        co_solvers.options['show_progress'] = False
        co_solvers.options['LPX_K_MSGLEV'] = 0
        co_solvers.options['MOSEK'] = {'mosek.iparam.log': 0}

        # the following options are for GLPK and are taken from http://www.maximal-usa.com/solvers/glpk.html
        # they should enhance performance for our large problems
        co_solvers.options['LPX_K_DUAL'] = 1 # use the dual simplex method with glpk because we likely have large problems
        co_solvers.options['LPX_K_PRESOL'] = 1 # use the linear presolver which was recommended to make things faster

        ncol = x.shape[1]
        xx = np.dot(x.transpose(), x)
        xy = np.dot(x.transpose(), y)
        c = co_matrix(np.ones(2 * ncol, dtype=float))
        G = np.vstack((-xx, xx))
        G = np.hstack((G, -G))
        G = np.vstack((G, -np.identity(2 * ncol)))
        G = co_matrix(G)
        h = co_matrix(np.hstack((sigma - xy, sigma + xy, np.zeros(2 * ncol))))

        # try each of the optimized solvers,
        # falling back on the unoptimized CVXOPT conelp solver when they aren't
        solvers = ('glpk', 'mosek', None)
        sol = None

        # keep the first working method around
        if cls.__i < 0:
            cls.__i = 0

        while sol is None and cls.__i < len(solvers):
            try:
                sol = co_solvers.lp(c, G.T, h, solver=solvers[cls.__i])
            except ValueError:
                cls.__i += 1
            if sol is not None and ('x' not in sol or sol['x'] is None):
                sol = None
                cls.__i += 1

        if sol is None or 'x' not in sol or sol['x'] is None:
            raise ValueError('problem has unbounded solution')

        sol = np.array([x for x in sol['x']])

        np.seterr(**nperr)

        # reset to defaults
        co_solvers.options.clear()

        return sol[:ncol] - sol[-ncol:]

    def learn(self, x, y):
        if not isinstance(x, np.ndarray):
            raise ValueError("x must be an numpy 2d array")

        if not isinstance(y, np.ndarray):
            raise ValueError("y must be an numpy 2d array")

        if x.ndim > 2:
            raise ValueError("x must be an 2d array")

        if x.shape[0] != y.shape[0]:
            raise ValueError("x and y are not aligned")

        # get default np error settings and set it to ignore divide errors
        nperr = np.geterr()
        np.seterr(divide='ignore')

        ncol = x.shape[1]

        if self.__m is None:
            self.__m = ncol

        sigma = self.__lam * np.sqrt(2.0 * np.log(ncol))
        mu = self.__glm.initialize(y)
        xbeta = self.__glm.link(mu)
        l = 0

        beta_old = np.empty(ncol, dtype=float)
        beta_old[:] = 2.0
        beta_new = np.empty(ncol, dtype=float)
        beta_new[:] = 1.0

        while sum(np.abs(beta_new - beta_old)) / (self.__tol + sum(np.abs(beta_old))) > self.__tol and l < self.__m:
            beta_old = np.copy(beta_new)
            v = np.mean(self.__glm.var(mu))
            z = xbeta + (y - mu) / v
            v = np.sqrt(v)
            zstar = z * v
            xstar = x * v
            # scale columns of Xstar to unit length
            colscale = np.array([np.sqrt(sum(np.square(x[:, j]))) for j in range(ncol)])
            xstar /= colscale
            beta_new = Dantzig.__linearoptim(xstar, zstar, sigma) / colscale
            xbeta = np.dot(x, beta_new)
            mu = self.__glm.invlink(xbeta)
            if sum(np.abs(beta_new)) < self.__tol:
                warn('beta all 0')
            l += 1

        # set anything below our tolerance to 0
        beta_new = beta_new * (np.abs(beta_new) > self.__tol)

        self.__L = l
        self.__beta = beta_new
        self.__selected = [i for i in range(ncol) if (self.__beta[i] != 0.0)]

        np.seterr(**nperr)

    def beta(self):
        return self.__beta

    def pred(self, x):
        if not isinstance(x, np.ndarray):
            raise ValueError("x must be an numpy 2d array")

        if x.ndim > 2:
            raise ValueError("x must be an 2d array")

        if x.shape[1] != self.__beta.shape[0]:
            print(x.shape[1], self.__beta.shape[0])
            raise ValueError("x and beta are not aligned")

        return np.dot(x, self.__beta)

    def selected(self):
        return self.__selected

    def steps(self):
        return self.__L
