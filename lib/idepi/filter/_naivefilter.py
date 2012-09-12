
from __future__ import division, print_function

from logging import getLogger
from operator import itemgetter
from time import clock

import numpy as np

from BioExt import Counter

from ..alphabet import Alphabet
from ._basefilter import BaseFilter
from ._filterutil import consgap_skip_columns
from ..logging import IDEPI_LOGGER


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
        self.__skip_cols = None
        self.__aln_meta = None
#         self.__run, self.__data, self.__colnames = False, None, None

    @staticmethod
    def __compute(alignment, alphabet, mincons, maxcons, maxgap,
                  refidx, refseq_offs, skip_func, skip_cols,
                  loop_defs):

        stream = posstream(alignment, alphabet, refidx)

        labels = []

        for i, pos in enumerate(stream):
            if (pos.maxcount / pos.total > maxcons or
            for char in (alphabet[j] for j in range(len(alphabet))):

            if pos.frac <=
        # PNGs binding site features: N-*[^PX-]-*[ST]

        colnames = BaseFilter._colnames(alignment, alphabet, refidx, refseq_offs, skip_dcols)
        data = NaiveFilter._filter(alignment, alphabet, skip_dcols, loop_defs, refidx)

        colnames.extend(k for k, _ in loop_defs)

        msg = "column labels don't equal number of columns: %d vs %d" % (len(colnames), data.shape[1])
        assert len(colnames) == data.shape[1], msg

        return colnames, data, aln_len, skip_dcols

    @staticmethod
    def _filter(alignment, alphabet, skip_cols, loop_defs, refidx=None):

        # if refidx exists, then we want to make sure
        # the reference isn't in the output
        nrow = len(alignment)
        if refidx is not None:
            nrow -= 1

        aln_len = alignment.get_alignment_length()
        stride = len(alphabet)
        ncol = aln_len * stride - len(skip_cols)
        data = np.zeros((nrow, ncol + len(loop_defs)), dtype=int)

        # build and save a list of offsets for each potential index
        offs = {}
        skipped = 0
        for i in range(aln_len * stride):
            if i in skip_cols:
                skipped += 1
            else:
                offs[i] = skipped

        b = clock()
        for j in range(aln_len):
            col = alignment[:, j]
            for i, char in enumerate(col):
                # handle skipping the reference
                r = i
                if refidx is not None:
                    if i == refidx:
                        continue
                    elif i > refidx:
                        r -= 1
                idx = j * stride + alphabet[char]
                if idx in skip_cols:
                    continue
                k = idx - offs[idx]
                data[r, k] = True

        k = ncol

        for start, end in (v for _, v in loop_defs):
            hvr = alignment[:, start:end]
            for i, row in enumerate(hvr):
                # handle skipping the reference
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

    def learn(self, alignment, refseq_offs, skip_cols=None):
        if skip_cols is None:
            skip_cols = consgap_skip_columns(
                alignment,
                self.__mincons,
                self.__maxcons,
                self.__maxgap,
                self.__refidx
            )
        colnames, data, self.__aln_len, self.__skip_cols = NaiveFilter.__compute(
            alignment,
            self.__alph,
            self.__mincons,
            self.__maxcons,
            self.__maxgap,
            self.__refidx,
            refseq_offs,
            self.__sfn,
            skip_cols,
            self.__loop_defs
        )
        return colnames, data

    def filter(self, alignment):
        if self.__aln_len is None or self.__skip_cols is None:
            raise RuntimeError('No NaiveFilter model computed')

        assert(alignment.get_alignment_length() == self.__aln_len)
        return NaiveFilter._filter(alignment, self.__alph, self.__skip_cols, self.__loop_defs)

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
