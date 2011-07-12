
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

    def __init__(self, alignment, alphabet=None, mincons=None, maxcons=None, maxgap=None, ref_id_func=None):

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

        self.__alph = alphabet
        self.__maxcons = maxcons
        self.__mincons = mincons
        self.__maxgap = maxgap
        self.__rfn = ref_id_func
#         self.__run, self.__data, self.__labels = False, None, None

    @staticmethod
    def __compute(alignment, alphabet, mincons, maxcons, maxgap, ref_id_func):

        seqtable = SeqTable(alignment, alphabet, ref_id_func)

        ignore_idxs = set()

        for i in xrange(seqtable.ncol):
            col = seqtable.cols[i]
            colsum = float(np.sum(col))
            if float(min(col)) / colsum > mincons or \
               float(max(col)) / colsum > maxcons or \
               float(col[alphabet['-']]) / colsum > maxgap:
                ignore_idxs.add(i)

        stride = len(alphabet)
        ncol = (len(seqtable) - len(ignore_idxs)) * stride 
        data = np.zeros((seqtable.nrow, ncol), dtype=bool)

        j = 0
        for i in xrange(seqtable.ncol):
            if i in ignore_idxs:
                continue
            data[:, j:(j+stride)] = seqtable.data[:, i]
            j += stride

        labels = BaseFilter._labels(alignment, alphabet, ref_id_func, ignore_idxs)

        return labels, data

    def filter(self, alignment):
        return NaiveFilter.__compute(
                alignment, self.__alph, self.__mincons, self.__maxcons, self.__maxgap, self.__rfn
            )

#     @property
#     def data(self):
#         if not self.__run:
#             raise RuntimeError('No naive filtering model computed')
#         return self.__data
# 
#     @property
#     def labels(self):
#         if not self.__run:
#             raise RuntimeError('No naive filtering model computed')
#         return self.__labels
