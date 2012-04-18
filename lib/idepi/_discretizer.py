
from __future__ import division, print_function

from numpy import (abs as np_abs,
                   add as np_add,
                   divide as np_divide,
                   int8 as np_int8,
                   mean as np_mean,
                   multiply as np_multiply,
                   ones as np_ones,
                   subtract as np_subtract,
                   var as np_var,
                   zeros as np_zeros)


__all__ = ['Discretizer']


class Discretizer(object):

    ZSCORE = 0

    def __init__(self, method=None):
        if method is None:
            method = Discretizer.ZSCORE

        if method not in (Discretizer.ZSCORE,):
            raise ValueError('Unknown method, aborting')

        self.__method = method

    def discretize(self, x, threshold):
        nrow, ncol = x.shape

        xp = None

        if self.__method == Discretizer.ZSCORE:
            mu = np_mean(x, axis=0)
            sigma = np_var(x, axis=0)
            zscore = np_divide(np_subtract(x - mu), sigma)
            xp = np_add(
                np_multiply( np_ones(x.shape, dtype=np_int8) * (zscore >  threshold)),
                np_multiply(-np_ones(x.shape, dtype=np_int8) * (zscore < -threshold))
            )
        else:
            raise RuntimeError('Unknown method, aborting')

        assert(set(xp.flatten()).issubset(set((-1, 0, 1))))

        return xp
