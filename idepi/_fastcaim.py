
import multiprocessing as mp
import numpy as np

# from _fakepool import FakePool


__all__ = ['FastCaim']


def _compute_caim(y, y_eye, intervals, nY, nI, tmp1=None):
    if tmp1 is None:
        tmp1 = np.zeros(y.shape, dtype=np.int32)
    res = 0.
    q, maxq = 0, 0
    for i in xrange(nI):
        tmp1[:] = intervals == i
        for j in xrange(nY):
            q = np.dot(tmp1, y_eye[j, :])
            if q > maxq:
                maxq = q
        M = np.sum(tmp1)
        res += np.nan_to_num(maxq / M)
    return res / nI


def _discretize(x, y):
    assert(x.shape == y.shape)

    nrow, = y.shape
    dtype = [('value', float), ('class', bool)]
    sortxy = np.zeros((nrow,), dtype=dtype)
    sortxy['value'] = x
    sortxy['class'] = y
    sortxy.sort(order='value')

    B = np.zeros((nrow - 1,), dtype=np.float64)
    nB = 0
    for i in xrange(1, nrow):
        # only insert boundaries if the class changes between two sorted variables
        a, b = sortxy[[i-1, i]]
        if a['class'] != b['class']:
            boundary = (a['value'] + b['value']) / 2.
            B[nB] = boundary
            nB += 1

    assert(nB < (nrow + 1))

    cidx = None
    caim = 0.
    innermaxcaim = 0.
    outermaxcaim = 0.
    k = 0 # use to count max iterations..

    D = np.zeros((nrow,), dtype=np.int32)
    intervals = np.zeros((nrow,), dtype=np.int32)

    included = set()

    nY = int(np.max(y)) + 1
    y_eye = np.zeros((nY, nrow), dtype=np.int32)
    for i in xrange(nY):
        y_eye[i, :] = y == i

    # make these here so that we don't constantly reinitialize arrays
    tmp1 = np.zeros((nrow,), dtype=np.int32)

    while True:
        for i in xrange(nB):
            # if i is in included, then we've already accounted for it in D, skip
            if i in included:
                continue
            # intervals is initialized to D
            intervals[:] = D
            # if x is greater than this proposed interval boundary,
            # increment its label -- this lets us add intervals
            # above and below the previous maximum boundary,
            # since compute_caim() is stateless.
            intervals += 1 * (x > B[i])
            caim = _compute_caim(y, y_eye, intervals, nY, k + 2, tmp1)
            if caim > innermaxcaim:
                innermaxcaim = caim
                cidx = i

        if innermaxcaim > outermaxcaim:
            outermaxcaim = innermaxcaim
            included.add(cidx)
            D += (x > B[cidx])
            k += 1
        else:
            break

    return D, outermaxcaim


class FastCaim(object):
    '''
    Implements the F-CAIM algorithm for binary classes from:
    Fast Class-Attribute Interdependence Maximization (CAIM) Discretization Algorithm
    by Kurgan and Cios, Oct 31 2010
    '''

    def __init__(self):
        pass

    def discretize(self, x, y):
        np_err = np.seterr(divide='ignore')

        nrow, ncol = x.shape
        self.__x = np.array(x.T, dtype=float, copy=True)
        self.__y = np.array(y, dtype=np.int32, copy=True).reshape(nrow)

        pool = mp.Pool(mp.cpu_count())

        res = [None] * ncol

        for j in xrange(ncol):
            res[j] = pool.apply_async(_discretize, (self.__x[j, :], self.__y))

        pool.close()
        pool.join()

        newx = np.zeros(self.__x.shape, dtype=int)
        caims = np.zeros((ncol,), dtype=float)

        for j in xrange(ncol):
            newx[j, :], caims[j] = res[j].get()

        np.seterr(**np_err)

        return newx.T
