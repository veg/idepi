
import numpy as np

from _alphabet import Alphabet
from _basefilter import BaseFilter
from _seqtable import SeqTable
from _util import is_HXB2


__all__ = ['NaiveFilter']


class NaiveFilter(BaseFilter):
    __DEFAULT_MIN_CONSERVATION = 1. # 100%
    __DEFAULT_MAX_CONSERVATION = 1. # 100%
    __DEFAULT_MAX_GAP_RATIO = 1.

    def __init__(self, alphabet=None, mincons=None, maxcons=None,
                 maxgap=None, ref_id_func=None, skip_func=None):

        if alphabet is None:
            alphabet = Alphabet()
        if maxcons is None:
            maxcons = NaiveFilter.__DEFAULT_MAX_CONSERVATION
        if mincons is None:
            mincons = NaiveFilter.__DEFAULT_MIN_CONSERVATION
        if maxgap is None:
            maxgap = NaiveFilter.__DEFAULT_MAX_GAP_RATIO
        if ref_id_func is None:
            ref_id_func = is_HXB2
        if skip_func is None:
            skip_func = lambda x: False

        self.__alph = alphabet
        self.__maxcons = maxcons
        self.__mincons = mincons
        self.__maxgap = maxgap
        self.__rfn = ref_id_func
        self.__sfn = skip_func
        self.__ignore_idxs = None
#         self.__run, self.__data, self.__colnames = False, None, None

    @staticmethod
    def __compute(alignment, alphabet, mincons, maxcons, maxgap,
                  ref_id_func, refseq_offs, skip_func):

        seqtable = SeqTable(alignment, alphabet, ref_id_func, skip_func)

        ignore_idxs = set()

        stride = len(alphabet)
        assert(len(seqtable.cols[0]) == stride)

        for i in xrange(seqtable.ncol):
            col = seqtable.cols[i]
            colsum = float(np.sum(col))
            idx = stride * i
            # kill perfectly conserved or empty columns 
            if colsum == 0. or \
               float(min(col)) / colsum > mincons or \
               float(max(col)) / colsum > maxcons or \
               float(max(col)) / colsum >= 100. or \
               float(col[alphabet['-']]) / colsum > maxgap:
                ignore_idxs.update(xrange(idx, idx + stride))
            else:
                ignore_idxs.update([idx + j for j in xrange(stride) if col[j] == 0.])

        colnames = BaseFilter._colnames(alignment, alphabet, ref_id_func, refseq_offs, ignore_idxs)
        data = NaiveFilter._filter(seqtable, ignore_idxs)

        return colnames, data, ignore_idxs

    @staticmethod
    def _filter(seqtable, ignore_idxs):
        stride = len(seqtable.cols[0])
        ncol = seqtable.ncol * stride - len(ignore_idxs)
        data = np.zeros((seqtable.nrow, ncol), dtype=bool)

        k = 0
        for i in xrange(seqtable.ncol):
            for j in xrange(stride):
                idx = i * stride + j
                if idx in ignore_idxs:
                    continue
                data[:, k] = seqtable.data[:, i, j]
                k += 1

        return data

    def learn(self, alignment, refseq_offs):
        colnames, data, self.__ignore_idxs = NaiveFilter.__compute(
            alignment, self.__alph, self.__mincons, self.__maxcons,
            self.__maxgap, self.__rfn, refseq_offs, self.__sfn
        )
        return colnames, data

    def filter(self, alignment):
        if not self.__ignore_idxs:
            raise RuntimeError('No NaiveFilter model computed')

        seqtable = SeqTable(alignment, self.__alph, self.__rfn, self.__sfn)
        assert(len(self.__alph) == len(seqtable.cols[0]))
        return NaiveFilter._filter(seqtable, self.__ignore_idxs)

#     @property
#     def data(self):
#         if not self.__run:
#             raise RuntimeError('No naive filtering model computed')
#         return self.__data
#
#     @property
#     def colnames(self):
#         if not self.__run:
#             raise RuntimeError('No naive filtering model computed')
#         return self.__colnames
