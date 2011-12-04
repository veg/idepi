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

from exceptions import AssertionError, NameError
from itertools import chain
from warnings import warn

from cvxopt import matrix as co_matrix
from cvxopt import solvers as co_solvers
import numpy as np
from mlpy import RidgeRegression, Lar, Lasso, ElasticNet, GradientDescent

# TODO: move these options into their proper functions (np.seterr will go everywhere...), co_solvers will go in Dantzig
# np.seterr(all='raise')
np.seterr(divide='ignore')

# silence the LP solvers
co_solvers.options['show_progress'] = False
co_solvers.options['LPX_K_MSGLEV'] = 0
co_solvers.options['MOSEK'] = {'mosek.iparam.log': 0}

# the following options are for GLPK and are taken from http://www.maximal-usa.com/solvers/glpk.html
# they should enhance performance for our large problems
co_solvers.options['LPX_K_DUAL'] = 1 # use the dual simplex method with glpk because we likely have large problems
co_solvers.options['LPX_K_PRESOL'] = 1 # use the linear presolver which was recommended to make things faster

__all__ = [
    'Glm',
    'LinearSvr',
    'Dantzig',
    'WrappedRegressor',
    'RidgeDantzig',
    'DoubleDantzig',
    'LarDantzig',
    'LassoDantzig',
    'RidgeLar',
    'RidgeLasso',
    'Regressor',
    'regressor_methods',
]

regressor_methods = {
    'ridge': 'RidgeRegression',
    'lar': 'Lar',
    'lasso': 'Lasso',
    'elasticnet': 'ElasticNet',
    'gd': 'GradientDescent',
    'dantzig': 'Dantzig',
#   'metadantzig': 'MetaDantzig', # can't define this here because
    'ridgedantzig': 'RidgeDantzig', # it's here but we don't really need it in light of MetaDantzig
    'doubledantzig': 'DoubleDantzig', # it's here but we don't really need it in light of MetaDantzig
    'lardantzig': 'LarDantzig',
    'lassodantzig': 'LassoDantzig',
    'ridgelar': 'RidgeLar',
    'ridgelasso': 'RidgeLasso',
    'lsvr': 'LinearSVR',
}


class Glm(object):

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


class LinearSvr(object):

    def __init__(self, tol=0.001, lam=1):
        pass

    def learn(self, x, y):
        pass

    def pred(self, x):
        pass

    def selected(self):
        pass

    def beta(self):
        pass


class Dantzig(object):
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


class WrappedRegressor(object):

    def __init__(self, method0, method1, **kwargs):

        methods = list(regressor_methods.keys())
        try:
            methods.pop(methods.index('metadantzig'))
        except ValueError:
            pass
        if method0 not in methods:
            raise ValueError('method `%s\' not in list of valid methods (%s)' % (method0, ', '.join(methods)))
        if method1 not in methods:
            raise ValueError('method `%s\' not in list of valid methods (%s)' % (method1, ', '.join(methods)))

        kwargs0 = {}
        kwargs1 = {}
        for k, v in kwargs.items():
            if k[-1] == '0':
                kwargs0[k[:-1]] = v
            else:
                if k[-1] == '1':
                    k = k[:-1]
                kwargs1[k] = v

        self.__first = eval(regressor_methods[method0])(**kwargs0)
        self.__second = eval(regressor_methods[method1])(**kwargs1)
        self.__beta0 = None
        self.__beta = None
        self.__idxs = None

    def learn(self, x, y):
        if not isinstance(x, np.ndarray):
            raise ValueError("x must be an numpy 2d array")

        if not isinstance(y, np.ndarray):
            raise ValueError("y must be an numpy 2d array")

        if x.ndim > 2:
            raise ValueError("x must be an 2d array")

        if x.shape[0] != y.shape[0]:
            raise ValueError("x and y are not aligned")

        self.__first.learn(x, y)
        self.__beta0 = self.__first.beta()

        self.__idxs = []
        for i in range(len(self.__beta0)):
            if self.__beta0[i] != 0.0:
                self.__idxs.append(i)
        x = x[:, self.__idxs]

        self.__second.learn(x, y)

        # reconstruct the full beta
        beta = np.zeros(len(self.__beta0))
        beta_ = self.__second.beta()
        for i in range(len(self.__idxs)):
            beta[self.__idxs[i]] = beta_[i]
        self.__beta = beta

    def beta0(self):
        return np.copy(self.__beta0)

    def beta(self):
        return np.copy(self.__beta)

    def pred(self, x):
        x = x[:, self.__idxs]
        return self.__second.pred(x)

    def selected(self):
        return self.__second.selected()

    def steps(self):
        if 'steps' in dir(self.__second):
            return self.__second.steps()
        warn('%s has no method `steps()\'' % self.__second.__class__.__name__)
        return None


class RidgeDantzig(WrappedRegressor):
    def __init__(self, tol=0.001, lam=1.0, alpha=0.0, **kwargs):
        if 'tol0' in kwargs:
            lam=kwargs['tol0']
            del kwargs['tol0']
        if 'lam0' in kwargs:
            lam=kwargs['lam0']
            del kwargs['lam0']
        super(RidgeDantzig, self).__init__(method0='dantzig', method1='ridge', tol0=tol, lam0=lam, alpha1=alpha, **kwargs)


class DoubleDantzig(WrappedRegressor):
    def __init__(self, tol0=0.001, lam0=1.0, tol=0.001, lam=0.001, **kwargs):
        super(DoubleDantzig, self).__init__(method0='dantzig', method1='dantzig', tol0=tol0, lam0=lam0, tol1=tol, lam1=lam, **kwargs)


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
        super(LarDantzig, self).__init__(method0='dantzig', method1='lar', tol0=tol, lam0=lam, m1=m, **kwargs)


class LassoDantzig(WrappedRegressor):
    def __init__(self, tol=0.001, lam=1.0, m=None, **kwargs):
        if m is None:
            raise ValueError('you must specify a maximum number of iterations m')
        if 'tol0' in kwargs:
            lam=kwargs['tol0']
            del kwargs['tol0']
        if 'lam0' in kwargs:
            lam=kwargs['lam0']
            del kwargs['lam0']
        super(LassoDantzig, self).__init__(method0='dantzig', method1='lasso', tol0=tol, lam0=lam, m1=m, **kwargs)


class RidgeLar(WrappedRegressor):
    def __init__(self, m, alpha=0.0, **kwargs):
        super(RidgeLar, self).__init__(method0='lar', method1='ridge', m0=m, alpha1=alpha, **kwargs)


class RidgeLasso(WrappedRegressor):
    def __init__(self, m, alpha=0.0, **kwargs):
        super(RidgeLasso, self).__init__(method='lasso', method1='ridge', m0=m, alpha1=alpha, **kwargs)


class Regressor(object):

    def __init__(self, data, method='ridgedantzig', *args, **kwargs):

        if method not in regressor_methods.keys():
            raise ValueError('method `%s\' not in list of valid methods (%s)' % (method, ', '.join(regressor_methods.keys())))

        self.__method = method
        self.__regressor = eval(regressor_methods[self.__method])(*args, **kwargs)
        self.__learned = False

        self.__feature_names = data.feature_names
        self.__x, self.__y = Regressor.__smldata_to_xy(data)
        self.__x_active, self.__y_active = self.__x.copy(), self.__y.copy()
        self.__x_inactive, self.__y_inactive = None, None

        # define the normalization constants
        # define these here to make pylint be quiet
        self.__xbar = None
        self.__xvar = None
        self.__ybar = None

        self.__folds = None
        self.__partition = None

    def __normalize(self):
        if self.__method not in ('gd',):
            __zero = pow(10.0, -9.0) # numerically close enough to 0.0
            ncol = self.__x_active.shape[1]

            # normalize y and validate
            self.__ybar = np.mean(self.__y_active)
            self.__y_active -= self.__ybar
            assert(np.abs(np.mean(self.__y_active)) < __zero)

            # normalize x and validate
            self.__xbar = np.zeros(ncol, dtype=float)
            self.__xvar = np.zeros(ncol, dtype=float)
            for j in range(ncol):
                self.__xbar[j] = np.mean(self.__x_active[:, j])
                self.__x_active[:, j] -= self.__xbar[j]
                assert(np.abs(np.mean(self.__x_active[:, j])) < __zero)
                if self.__method not in ('ridge', 'gd'):
                    self.__xvar[j] = np.sqrt(sum(pow(self.__x_active[:, j], 2.0)))
                    if self.__xvar[j] != 0.0:
                        self.__x_active[:, j] /= self.__xvar[j]
                        try:
                            assert(np.abs(sum([pow(i, 2.0) for i in self.__x_active[:, j]]) - 1.0) < __zero)
                        except AssertionError as e:
                            print('\u03c3: %.4g, \u03a3x\u00b2: %.4g' % (self.__xvar[j], sum([pow(i, 2.0) for i in self.__x_active[:, j]])))

    def partition(self, l, folds):
        npf = int(np.floor(l / folds)) # num per fold
        r = len % folds
        p = list(chain(*([i] * npf for i in range(folds)))) + list(range(r))
        np.random.shuffle(p)
        assert(len(p) == self.__x.shape[0])
        self.__folds = folds
        self.__partition = p

    def mask(self, fold):
        assert(0 <= fold and fold < self.__folds)
        active = [i for i in range(self.__x.shape[0]) if self.__partition[i] != fold]
        inactive = [i for i in range(self.__x.shape[0]) if self.__partition[i] == fold]
        # don't need to issue a copy() here
        self.__x_active = self.__x[active, :]
        self.__y_active = self.__y[active]
        self.__x_inactive = self.__x[inactive, :]
        self.__y_inactive = self.__y[inactive]

    def unmask(self):
        self.__x_active, self.__y_active = self.__x.copy(), self.__y.copy()
        self.__x_inactive, self.__y_inactive = None, None

    @staticmethod
    def __smldata_to_xy(data):
        ncol = len(data.feature_names)

        x = []
        y = []

        for row in data:
            x.append([row.features[j] if j in row.features else 0. for j in range(ncol)])
            y.append(row.value)

        x = np.array(x).astype(float)
        x = np.hstack((np.ones((x.shape[0], 1), dtype=float), x)) # add a constant term (y-intercept) as the 0th column
        y = np.array(y).astype(float)

        return x, y

    def __normalize_for_predict(self, x, y):
        if self.__xbar is None or self.__xvar is None or self.__ybar is None:
            raise NameError('normalization constants (xbar, xvar, yvar) are unset')

        ncol = x.shape[1]

        if self.__method not in ('gd',):
            y -= self.__ybar # LAR, LASSO, ElasticNet: responses have mean 0
            # skip index 0 because it's our y-intercept
            for j in range(1, ncol):
                x[:, j] -= self.__xbar[j] # LAR, LASSO, ElasticNet: covariates have mean 0
                if self.__method not in ('ridge', 'gd'):
                    if self.__xvar[j] != 0.0:
                        x[:, j] /= self.__xvar[j] # LAR, LASSO: covariates have unit length

        return x, y

    @property
    def intercept(self):
        if not self.__learned:
            Regressor.learn(self)
        return self.__regressor.beta()[0]

    @property
    def weights(self):
        if not self.__learned:
            Regressor.learn(self)
        beta = self.__regressor.beta()
        # ignore the y-intercept
        return dict(zip(self.__feature_names, beta[1:]))

    @property
    def selected(self):
        if not self.__learned:
            Regressor.learn(self)
        selected = self.__regressor.selected()
        # -1 refers to the y-intercept
        return np.array([s - 1 for s in selected])

    def learn(self):
        self.__regressor.learn(self.__x_active, self.__y_active)
        self.__learned = True

    def predict(self, data):
        x, y = Regressor.__smldata_to_xy(data)
        x, y = Regressor.__normalize_for_predict(self, x, y)
        if not self.__learned:
            Regressor.learn(self)
        return self.__regressor.pred(x)

    def test(self, data):
        if not self.__learned:
            Regressor.learn(self)
        x, y = Regressor.__smldata_to_xy(data)
        x, y = Regressor.__normalize_for_predict(self, x, y)
        yhat = self.__regressor.pred(x)
        sse = sum(pow(y - yhat, 2.0))
        ybar = np.mean(y)
        sst = sum(pow(y - ybar, 2.0))
        r2 = 1.0 - (sse / sst)
        nless1 = len(y) - 1
        p = len(1 for i in self.weights.values() if i != 0.0) - 1 # - 1 to avoid counting the constant term
        mse = sse / (nless1 - p) # count the the full N
        rmse = np.sqrt(mse)
        rbar2 = 1.0 - (1.0 - r2) * nless1 / (nless1 - p)
        return {
            'R\u0304\u00b2   ': rbar2,
            'R\u00b2   ': r2,
            # u'SSE  ': sse,
            # u'MSE  ': mse,
            'RMSE ': rmse,
        }
