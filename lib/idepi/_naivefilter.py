
from collections import Counter

import numpy as np

from ._alphabet import Alphabet
from ._basefilter import BaseFilter
from ._util import is_HXB2


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
        self.__aln_len = None
        self.__ign_idxs = None
#         self.__run, self.__data, self.__colnames = False, None, None

    @staticmethod
    def __compute(alignment, alphabet, mincons, maxcons, maxgap,
                  ref_id_func, refseq_offs, skip_func):

        aln_len = alignment.get_alignment_length()
        ign_idxs = set()
        stride = len(alphabet)
        for i in range(aln_len):
            counts = Counter(alignment[:, i])
            col = [v for k, v in counts if k != Alphabet.SPACE]
            colsum = sum(col)
            gapdiv = colsum + (counts[Alphabet.SPACE] if Alphabet.SPACE in counts else 0)
            idx = stride * i
            # kill perfectly conserved or empty columns 
            if (colsum == 0. or
                min(col) / colsum > mincons or
                max(col) / colsum > maxcons or
                max(col) / colsum >= 1. or
                (Alphabet.SPACE in counts and counts[Alphabet.SPACE] / gapdiv > maxgap)):
                ign_idxs.update(range(idx, idx + stride))
            else:
                ign_idxs.update(idx + j for j in range(stride) if col[j] == 0.)

        colnames = BaseFilter._colnames(alignment, alphabet, ref_id_func, refseq_offs, ign_idxs)
        data = NaiveFilter._filter(alignment, alphabet, ign_idxs)

        assert(len(colnames) == data.shape[1])

        return colnames, data, aln_len, ign_idxs

    @staticmethod
    def _filter(alignment, alphabet, ign_idxs):
        aln_len = alignment.get_alignment_length()
        stride = len(alphabet)
        ncol = aln_len * stride - len(ign_idxs)
        data = np.zeros((len(alignment), ncol), dtype=bool)

        k = 0
        for i in range(aln_len):
            col = alignment[:, i]
            for j, char in enumerate(alphabet):
                idx = i * stride + j
                if idx in ign_idxs:
                    continue
                data[:, k] = (p == char for p in col)
                k += 1

        return data

    def learn(self, alignment, refseq_offs):
        colnames, data, self.__aln_len, self.__ign_idxs = NaiveFilter.__compute(
            alignment, self.__alph, self.__mincons, self.__maxcons,
            self.__maxgap, self.__rfn, refseq_offs, self.__sfn
        )
        return colnames, data

    def filter(self, alignment):
        if self.__aln_len is None or self.__ign_idxs is None:
            raise RuntimeError('No NaiveFilter model computed')

        assert(alignment.get_alignment_length() == self.__aln_len)
        return NaiveFilter._filter(alignment, self.__alph, self.__ign_idxs)

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
