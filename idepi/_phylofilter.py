
import json
from math import floor
from os import close, getpid, remove
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

    def __init__(self, alphabet=None, batchfile=None, ref_id_func=None, skip_func=None):
        if batchfile is None:
            batchfile = join(dirname(realpath(__file__)), '..', 'res', 'hyphy', 'CorrectForPhylogeny.bf')

        if not exists(batchfile):
            raise ValueError('Please pass a valid (and existing) batchfile to PhyloFilter()')

        if alphabet is None:
            alphabet = Alphabet()
        if ref_id_func is None:
            ref_id_func = is_HXB2
        if skip_func is None:
            skip_func = lambda x: False

        fd, self.__inputfile = mkstemp(); close(fd)

        self.__alph = alphabet
        self.__batchfile = batchfile
        self.__rfn = ref_id_func
        self.__sfn = skip_func
#         self.__run, self.__data, self.__colnames = False, None, None

    def __del__(self):
        for file in (self.__inputfile,):
            if file and exists(file):
                remove(file)

    @staticmethod
    def __compute(alignment, alphabet, batchfile, inputfile, ref_id_func, refseq_offs, skip_func, hyphy=None):
        if hyphy is None:
            hyphy = HyPhy()

        refseq = None
        for row in alignment:
            r = apply(ref_id_func, (row.id,))
            if r and refseq is None:
                refseq = row
            elif r:
                raise RuntimeError('Reference sequence found twice!?!?!?!')

        alignment = [row for row in alignment if \
                     not apply(ref_id_func, (row.id,)) and \
                     not apply(skip_func, (row.id,))
                    ]

        with open(inputfile, 'w') as fh:
            SeqIO.write(alignment, fh, 'fasta')

        HyPhy.execute(hyphy, batchfile, (inputfile,))

        order = HyPhy.retrieve(hyphy, 'order', HyPhy.STRING).strip(',').split(',')
#         _ids = HyPhy.retrieve(hyphy, 'ids', HyPhy.MATRIX)
        msm = HyPhy.retrieve(hyphy, 'marginalSupportMatrix', HyPhy.MATRIX)
        bim = HyPhy.retrieve(hyphy, 'binaryIdentitiesMatrix', HyPhy.MATRIX)
        nspecies = int(HyPhy.retrieve(hyphy, 'numSpecies', HyPhy.NUMBER))
        nsites = int(HyPhy.retrieve(hyphy, 'numSites', HyPhy.NUMBER))
        nchars = int(HyPhy.retrieve(hyphy, 'numChars', HyPhy.NUMBER))

        assert(len(order) == nchars)
        assert(msm.mCols == (nsites * nchars))
        assert(msm.mRows == nspecies)
        assert(bim.mCols == (nsites * nchars))
        assert(bim.mRows == nspecies)
#         assert(_ids.mRows == 0)
#         ids = [_ids.MatrixCell(0, i) for i in xrange(_ids.mCols)]

        ncol = nsites * len(alphabet)
        custom_type = np.dtype([('b', bool), ('p', float)]) # 'b' for binary identity, 'p' for probability | phylogeny
        tmp = np.zeros((nspecies, ncol), dtype=custom_type) #, order='F') # use column-wise order in memory

        # cache the result for each stride's indexing into the alphabet
        alphidx = [alphabet[order[i]] for i in xrange(len(order))]

        for i in xrange(nspecies):
            for j in xrange(nsites):
                for k in xrange(nchars):
                    # we map j from HyPhy column order into self.__alph column order
                    # by getting at the MSA column (j / len(order)), multiplying by
                    # the self.__alph stride (len(self.__alph)), and then finally adding
                    # the alphabet-specific index (alphidx[r])
                    l = (j * len(alphabet)) + alphidx[k]
                    # 'b' for binary identity and 'p' for probability | phylogeny
                    tmp[i, l] = (bim.MatrixCell(i, j*nchars + k), msm.MatrixCell(i, j*nchars + k))

        # np.save('phylofilt.%d' % getpid(), tmp)

        colsum = np.sum(tmp[:, :]['b'], axis=0)
        idxs = [i for i in xrange(ncol) if colsum[i] != 0]
        ignore_idxs = set([i for i in xrange(ncol) if colsum[i] == 0])

        data = tmp[:, idxs]

        np.savez('phylo.npz', {'data': data})

#         data = np.zeros((mat.mRows, ncol - len(ignore_idxs)), dtype=float)
#
#         j = 0
#         for i in xrange(ncol):
#             if i in ignore_idxs:
#                 continue
#             data[:, j] = tmp[:, i]
#             j += 1

        if refseq is not None:
            alignment.append(refseq)

        colnames = BaseFilter._colnames(alignment, alphabet, ref_id_func, refseq_offs, ignore_idxs)

        with open('phylo.json', 'w') as fh:
            import json
            json.dump(colnames, fh)

        # make sure that the columns do line up
        assert(len(colnames) == data.shape[1])

        # return ids, mat, order, colnames
        return colnames, data

    @staticmethod
    def filter(alignment):
        raise RuntimeError('PhyloFilter does not yet support the learn() and filter() model of its NaiveFilter brother')

    def learn(self, alignment, refseq_offs):
        return PhyloFilter.__compute(
            alignment, self.__alph, self.__batchfile, self.__inputfile, self.__rfn, refseq_offs, self.__sfn
        )

#     @property
#     def data(self):
#         if not self.__run:
#             raise RuntimeError('No phylofiltering model computed')
#         return self.__data
#
#     @property
#     def colnames(self):
#         if not self.__run:
#             raise RuntimeError('No phylofiltering model computed')
#         return self.__colnames
