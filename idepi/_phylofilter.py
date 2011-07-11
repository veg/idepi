
import json
from os import close, remove
from os.path import dirname, exists, join, realpath
from tempfile import mkstemp

from numpy import zeros

from Bio import SeqIO

from _alphabet import Alphabet
from _hyphy import HyPhy


__all__ = ['PhyloFilter']


class PhyloFilter(object):

    def __init__(self, seqrecords, alphabet=None, batchfile=None):
        if batchfile is None:
            batchfile = join(dirname(realpath(__file__)), '..', 'res', 'CorrectForPhylogeny.bf')

        if not exists(batchfile):
            raise ValueError('Please pass a valid (and existing) batchfile to PhyloFilter()')

        if alphabet is None:
            alphabet = Alphabet()

        fd, self.__inputfile = mkstemp(); close(fd)

        self.__seqrecords = seqrecords
        self.__alph = alphabet
        self.__batchfile = batchfile
        self.__commands = commands
        self.__hyphy = HyPhy()

        self.__ids, self.__mat, self.__ord = PhyloFilter.__run(self, seqrecords)

        self.__run = True

    def __del__(self):
        for file in (self.__inputfile,):
            if file and exists(file):
                remove(file)

    def __get_value(self, variable, type):
        _res = self.__hyphy.AskFor(variable)
        if type not in (HyPhy.MATRIX, HyPhy.NUMBER, HyPhy.STRING):
            raise ValueError('Unknown type supplied: please use one of PhyloFilter.{MATRIX,NUMBER,STRING}')
        if (self.__hyphy.CanICast(_res, type)):
            res = self.__hyphy.CastResult(_res, type)
            if type == HyPhy.STRING:
                return res.castToString().sData
            elif type == HyPhy.NUMBER:
                return res.castToNumber().nValue
            elif type == HyPhy.MATRIX:
                return res.castToMatrix()
            else:
                # dead code, we assume
                assert(0)
        else:
            raise RuntimeError('Cast failed in HyPhy, assume an incorrect type was supplied for variable `%s\'' % variable)

    def __run(self, seqrecords):
        with open(self.__inputfile, 'w') as fh:
            SeqIO.write(seqrecords, fh, 'fasta')

        self.__hyphy.ExecuteBF('ExecuteAFile("%s", { "0": "%s" })' % (self.__batchfile, self.__inputfile))

        _ids = PhyloFilter.__get_value(self, 'ids', HyPhy.MATRIX)
        _mat = PhyloFilter.__get_value(self, 'data', HyPhy.MATRIX)
        order = PhyloFilter.__get_value(self, 'order', HyPhy.STRING).split(',')

        assert(_ids.mRows == 0)

        ids = [_ids.MatrixCell(0, i) for i in xrange(_ids.mCols)]
        mat = zeros((_mat.mRows, _mat.mCols), dtype=float)

        for i in xrange(_mat.mRows):
            for j in xrange(_mat.mCols):
                mat[i, j] = _mat.MatrixCell(i, j)

        return ids, mat, order

    def names(self, ref_id_func):
        if not self.__run:
            raise RuntimeError('No phylofiltering model computed')

        ref = None
        for r in self.__seqrecords:
            if apply(ref_id_func, (r.id,)):
                ref = str(r.seq)

        if ref is None:
            raise RuntimeError('No reference sequence found, aborting')

        for i in xrange(len(ref)):
            pass
