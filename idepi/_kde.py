
from collections import namedtuple
from math import ceil, floor, pi, sqrt

import numpy as np


__all__ = ['Interval', 'KernelDensityEstimator']


Interval = namedtuple('Interval', ['min', 'max'])


def KernelDensityEstimator(data, n=2 ** 14, interval=None):
    '''
    Implements the kernel density estimator described in:
    Kernel density estimation via diffusion
    Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
    Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
    '''

    def bisect(f, a, b, args):
        tol = 4.4408920985006262e-16 
        fa = apply(f, (a,)+args)
        fb = apply(f, (b,)+args)
        if cmp(fa, 0) == cmp(fb, 0):
            raise RuntimeError('f(a) and f(b) must have different signs')
        while abs(fa) > tol and abs(fb) > tol:
            ab = (a + b) / 2.
            fab = apply(f, (ab,)+args)
            if cmp(fa, 0) != cmp(fab, 0):
                b = ab
                fb = fab
            elif cmp(fb, 0) != cmp(fab, 0):
                a = ab
                fa = fab
            else:
                raise RuntimeError('Something fishy happened during bisection')
        if abs(fa) < tol and abs(fb) < tol:
            return a if min(abs(fa), abs(fb)) == abs(fa) else b
        elif abs(fa) < tol:
            return a
        elif abs(fb) < tol:
            return b
        else:
            raise RuntimeError('Something fishy happened during bisection')

    def unique(lst):
        seen = set()
        return [x for x in lst if x not in seen and not seen.add(x)]

    def dct1d(data):
        nrow, = data.shape
        weights = 2 * np.exp(-1j * np.arange(nrow) * np.pi / (2. * nrow)) 
        weights[0] -= 1.
        data2 = data[range(0, nrow, 2) + range(nrow-1, 0, -2)]
        return np.real(np.multiply(weights, np.fft.fft(data2)))

    def idct1d(data):
        nrow, = data.shape
        weights = nrow * np.exp(1j * np.arange(nrow) * np.pi / (2. * nrow))
        data2 = np.real(np.fft.ifft(np.multiply(weights, data)))
        out = np.zeros((nrow,), dtype=float)
        out[range(0, nrow, 2)]    = data2[:(nrow / 2)]
        out[range(nrow-1, 0, -2)] = data2[(nrow / 2) + 1:]
        return out

    def compf(t, s, I, a2):
        return 2. * (pi ** (2. * s)) * np.sum(np.multiply(np.multiply(I ** s, a2), np.exp(-I * (np.pi ** 2) * t)))
    
    def fixed_point(t, N, I, a2):
        l = 7
        f = compf(t, l, I, a2)
        for s in xrange(l-1, 1, -1):
            K0 = np.prod(np.arange(1, 2 * s, 2)) / sqrt(2. * pi)
            const = (1. + (0.5 ** (s + 0.5))) / 3.
            time = (2. * const * K0 / N / f) ** (2. / (3 + (2 * s)))
            f = compf(time, s, I, a2)
        return t - ((2. * N * sqrt(pi) * f) ** (-2. / 5))

    if interval is None:
        dmax = max(data)
        dmin = min(data)
        r = dmax - dmin
        interval = Interval(min=dmin - (r / 10.), max=dmax + (r / 10.))

    if type(interval) != Interval:
        raise ValueError('interval must be of type Interval')

    nperr = np.seterr(under='ignore')

    n = 2 ** ceil(np.log2(n))

    R = interval.max - interval.min
    dx = 1. * R / (n - 1)
    N = len(unique(data))

    hist, mesh = np.histogram(data, n - 1, range=(interval.min, interval.max))
    initial_data = np.zeros(mesh.shape, dtype=float)
    initial_data[:-1] = hist
    initial_data /= N
    initial_data = initial_data / np.sum(initial_data)

    a = dct1d(initial_data)

    I = np.arange(1, n) ** 2
    a2 = (a[1:] / 2) ** 2

    if len(I) != len(a2):
        raise RuntimeError('Lengths of `I\' and `a2\' are different, respectively: %d, %d' % (len(I), len(a2)))

    try:
        # raise ImportError()
        from scipy.optimize import brentq
    except ImportError:
        brentq = bisect

    try:
        t_star = brentq(fixed_point, 0., 0.1, args=(N, I, a2))
    except ValueError:
        t_star = 0.28 * (N ** (-2. / 5))

    bandwidth = sqrt(t_star) * R

    a_t = np.multiply(a, np.exp(-(np.arange(n) ** 2) * (np.pi ** 2) * t_star / 2))

    density = idct1d(a_t) / R

    f = compf(t_star, 1, I, a2)
    t_cdf = (sqrt(pi) * f * N) ** (-2. / 3)
    a_cdf = np.multiply(a, np.exp(-(np.arange(n) ** 2) * (np.pi ** 2) * t_cdf / 2))
    pdf = idct1d(a_cdf) * (dx / R)
    cdf = np.cumsum(pdf)
    bandwidth_cdf = sqrt(t_cdf) * R
    
    np.seterr(**nperr)

    return (bandwidth, density, mesh, pdf, cdf)


def main():
    # from _data import DATA as data
    from time import time

    d1 = np.random.randn(100) + 5
    d2 = np.random.randn(100) * 2 + 35
    d3 = np.random.randn(100) + 55

    data = np.concatenate((d1, d2, d3))

    begin = time()

    bandwidth, density, mesh, pdf, cdf = KernelDensityEstimator(data, 2 ** 14) 

    runtime = time() - begin

    print runtime, bandwidth, len(density), np.sum(density), len(mesh), sum(pdf), cdf[-1]

    import matplotlib.pyplot as plt

    # plt.plot(mesh, density)
    plt.plot(mesh, pdf)
    plt.plot(mesh, cdf / 1000.)
    plt.show()

    return 0


if __name__ == '__main__':
    from sys import exit
    exit(main())
