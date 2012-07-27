
from __future__ import division, print_function

from logging import getLogger
from operator import itemgetter
from time import clock

import numpy as np

from BioExt import Counter

from ._alphabet import Alphabet
from ._basefilter import BaseFilter
from ._logging import IDEPI_LOGGER


__all__ = ['NaiveFilter']


class NaiveFilter(BaseFilter):
    __DEFAULT_MIN_CONSERVATION = 1. # 100%
    __DEFAULT_MAX_CONSERVATION = 1. # 100%
    __DEFAULT_MAX_GAP_RATIO = 1.

    def __init__(self, alphabet=None, mincons=None, maxcons=None,
                 maxgap=None, refidx=None, skip_func=None,
                 loop_defs=None):

        if alphabet is None:
            alphabet = Alphabet()
        if maxcons is None:
            maxcons = NaiveFilter.__DEFAULT_MAX_CONSERVATION
        if mincons is None:
            mincons = NaiveFilter.__DEFAULT_MIN_CONSERVATION
        if maxgap is None:
            maxgap = NaiveFilter.__DEFAULT_MAX_GAP_RATIO
        if loop_defs is None:
            loop_defs = []
        if skip_func is None:
            skip_func = lambda x: False

        self.__alph = alphabet
        self.__maxcons = maxcons
        self.__mincons = mincons
        self.__maxgap = maxgap
        self.__refidx = refidx
        self.__sfn = skip_func
        self.__loop_defs = loop_defs
        self.__aln_len = None
        self.__ign_idxs = None
#         self.__run, self.__data, self.__colnames = False, None, None

    @staticmethod
    def __compute(alignment, alphabet, mincons, maxcons, maxgap,
                  refidx, refseq_offs, skip_func,
                  loop_defs):

        b = clock()
        aln_len = alignment.get_alignment_length()
        ign_idxs = set()
        stride = len(alphabet)
        idx_chars = sorted(((v, k) for k, v in alphabet.todict().items()), key=itemgetter(0))
        for i in range(aln_len):
            idx_i = stride * i
            if idx_i in ign_idxs:
                continue
            for j in range(aln_len):
                counts_j = Counter(alignment[k, j] for k in range(len(alignment)) if k != refidx)
                col = [v for k, v in counts_j.items() if k != Alphabet.SPACE]
                colsum_nogaps = sum(col)
                colsum_withgaps = colsum_nogaps + (counts_j[Alphabet.SPACE] if Alphabet.SPACE in counts_j else 0)
                idx_j = stride * j
                # kill perfectly conserved or empty columns
                if (colsum_nogaps == 0. or
                    min(col) / colsum_nogaps > mincons or
                    max(col) / colsum_nogaps > maxcons or
                    max(col) / colsum_nogaps >= 1. or
                    (Alphabet.SPACE in counts_j and counts_j[Alphabet.SPACE] / colsum_withgaps > maxgap)):
                    ign_idxs.update(range(idx_j, idx_j + stride))
                else:
                    ign_idxs.update(idx_j + k for k, char in idx_chars if char not in counts_j or counts_j[char] == 0)

        getLogger(IDEPI_LOGGER).debug('finished learning a filter, took %.3f' % (clock() - b))

        # PNGs binding site features: N-*[^PX-]-*[ST]

        colnames = BaseFilter._colnames(alignment, alphabet, refidx, refseq_offs, ign_idxs)
        data = NaiveFilter._filter(alignment, alphabet, ign_idxs, loop_defs, refidx)

        colnames.extend(k for k, _ in loop_defs)

        msg = "column labels don't equal number of columns: %d vs %d" % (len(colnames), data.shape[1])
        assert len(colnames) == data.shape[1], msg

        return colnames, data, aln_len, ign_idxs

    @staticmethod
    def _filter(alignment, alphabet, ign_idxs, loop_defs, refidx=None):

        nrow = len(alignment)
        if refidx is not None:
            nrow -= 1

        aln_len = alignment.get_alignment_length()
        stride = len(alphabet)
        ncol = aln_len * stride - len(ign_idxs)
        data = np.zeros((nrow, ncol + len(loop_defs)), dtype=int)

        # build and save a list of offsets for each potential index
        offs = {}
        skipped = 0
        for i in range(aln_len * stride):
            if i in ign_idxs:
                skipped += 1
            else:
                offs[i] = skipped

        b = clock()
        for j in range(aln_len):
            col = alignment[:, j]
            for i, char in enumerate(col):
                r = i
                if refidx is not None:
                    if i == refidx:
                        continue
                    elif i > refidx:
                        r -= 1
                idx = j * stride + alphabet[char]
                if idx in ign_idxs:
                    continue
                k = idx - offs[idx]
                data[r, k] = True

        k = ncol

        for start, end in (v for _, v in loop_defs):
            hvr = alignment[:, start:end]
            for i, row in enumerate(hvr):
                r = i
                if refidx is not None:
                    if i == refidx:
                        continue
                    elif i > refidx:
                        r -= 1
                data[r, k] = sum(1 for p in row if p != Alphabet.SPACE)
            k += 1

        getLogger(IDEPI_LOGGER).debug('finished building a data matrix, took %.3fs' % (clock() - b))

        return data

    def learn(self, alignment, refseq_offs):
        colnames, data, self.__aln_len, self.__ign_idxs = NaiveFilter.__compute(
            alignment, self.__alph, self.__mincons, self.__maxcons,
            self.__maxgap, self.__rfn, refseq_offs, self.__sfn,
            self.__loop_defs
        )
        return colnames, data

    def filter(self, alignment):
        if self.__aln_len is None or self.__ign_idxs is None:
            raise RuntimeError('No NaiveFilter model computed')

        assert(alignment.get_alignment_length() == self.__aln_len)
        return NaiveFilter._filter(alignment, self.__alph, self.__ign_idxs, self.__loop_defs)

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
