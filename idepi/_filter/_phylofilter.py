
import json
from os import close, remove
from os.path import dirname, exists, join, realpath
from tempfile import mkstemp

from np import zeros

from Bio import SeqIO

from _hyphy import HyPhy


__all__ = ['PhyloFilter']


class PhyloFilter(object):
    MATRIX = HyPhy.THYPHY_TYPE_MATRIX
    NUMBER = HyPhy.THYPHY_TYPE_NUMBER
    STRING = HyPhy.THYPHY_TYPE_STRING

    def __init__(self, seqrecords, batchfile=None):
        if batchfile is None:
            batchfile = join(dirname(realpath(__file__)), '..', 'res', 'CorrectForPhylogeny.bf')

        if not exists(batchfile):
            raise ValueError('Please pass a valid (and existing) batchfile to PhyloFilter()')

        fd, self.__inputfile = mkstemp(); close(fd)

        self.__seqrecords = seqrecords
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
        if type not in (PhyloFilter.MATRIX, PhyloFilter.NUMBER, PhyloFilter.STRING):
            raise ValueError('Unknown type supplied: please use one of PhyloFilter.{MATRIX,NUMBER,STRING}')
        if (self.__hyphy.CanICast(_res, type):
            res = self.__hyphy.CastResult(_res, type)
            if type == PhyloFilter.STRING:
                return res.castToString().sData
            elif type == PhyloFilter.NUMBER:
                return res.castToNumber().nValue
            elif type == PhyloFilter.MATRIX:
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

        _ids = PhyloFilter.__get_value(self, 'ids', PhyloFilter.MATRIX)
        _mat = PhyloFilter.__get_value(self, 'data', PhyloFilter.MATRIX)
        order = PhyloFilter.__get_value(self, 'order', PhyloFilter.STRING).split(',')

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
