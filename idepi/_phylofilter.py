
import json
from math import floor
from os import close, remove
from os.path import dirname, exists, join, realpath
from tempfile import mkstemp

import numpy as np

from Bio import SeqIO

from _alphabet import Alphabet
from _basefilter import BaseFilter
from _hyphy import HyPhy
from _util import is_HXB2


__all__ = ['PhyloFilter']


class PhyloFilter(BaseFilter):

    def __init__(self, alphabet=None, batchfile=None, ref_id_func=None):
        if batchfile is None:
            batchfile = join(dirname(realpath(__file__)), '..', 'res', 'CorrectForPhylogeny.bf')

        if not exists(batchfile):
            raise ValueError('Please pass a valid (and existing) batchfile to PhyloFilter()')

        if alphabet is None:
            alphabet = Alphabet()
        if ref_id_func is None:
            ref_id_func = is_HXB2

        fd, self.__inputfile = mkstemp(); close(fd)

        self.__alph = alphabet
        self.__batchfile = batchfile
        self.__rfn = ref_id_func
#         self.__run, self.__data, self.__labels = False, None, None

    def __del__(self):
        for file in (self.__inputfile,):
            if file and exists(file):
                remove(file)

    @staticmethod
    def __compute(alignment, alphabet, batchfile, inputfile, ref_id_func, hyphy=None):
        if hyphy is None:
            hyphy = HyPhy()

        with open(inputfile, 'w') as fh:
            SeqIO.write(alignment, fh, 'fasta')

        HyPhy.execute(hyphy, batchfile, (inputfile,))

        _ids  = HyPhy.retrieve(hyphy, 'ids', HyPhy.MATRIX)
        _mat  = HyPhy.retrieve(hyphy, 'data', HyPhy.MATRIX)
        order = HyPhy.retrieve(hyphy, 'order', HyPhy.STRING).strip(',').split(',')

        assert(_ids.mRows == 0)

        ids = [_ids.MatrixCell(0, i) for i in xrange(_ids.mCols)]

        ncol = _mat.mCols / len(order) * len(alphabet)
        mat = np.zeros((_mat.mRows, ncol), dtype=float, order='F') # use column-wise order in memory

        # cache the result for each stride's indexing into the alphabet
        alphidx = []
        for i in xrange(len(order)):
            alphidx.append(alphabet[order[i]])

        for j in xrange(_mat.mCols):
            # we map j from HyPhy column order into self.__alph column order
            # by getting at the MSA column (j / len(order)), multiplying by
            # the self.__alph stride (len(self.__alph)), and then finally adding
            # the alphabet-specific index (alphidx[r])
            q = int(floor(j / len(order))) # quotient
            r = j % len(order) # remainder
            k = (q * len(alphabet)) + alphidx[r]
            for i in xrange(_mat.mRows):
                mat[i, k] += _mat.MatrixCell(i, j)

        ignore_idxs = set()
        colsum = np.sum(mat, axis=0)
        for j in xrange(ncol):
            if colsum[j] == 0.:
                ignore_idxs.add(j)

        labels = BaseFilter._labels(alignment, alphabet, ref_id_func, ignore_idxs)

        # return ids, mat, order, labels
        return labels, mat

    def filter(self, alignment):
        return PhyloFilter.__compute(
                alignment, self.__alph, self.__batchfile, self.__inputfile, self.__rfn
            )

#     @property
#     def data(self):
#         if not self.__run:
#             raise RuntimeError('No phylofiltering model computed')
#         return self.__data
# 
#     @property
#     def labels(self):
#         if not self.__run:
#             raise RuntimeError('No phylofiltering model computed')
#         return self.__labels
