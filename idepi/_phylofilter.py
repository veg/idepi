
import json
from os import close, remove
from os.path import dirname, exists, join, realpath
from tempfile import mkstemp

from np import zeros

from Bio import SeqIO

from _hyphy import HyPhy


__all__ = ['PhyloFilter']


class PhyloFilter(object):

    def __init__(self, batchfile=None):
        if batchfile is None:
            batchfile = join(dirname(realpath(__file__)), '..', 'res', 'CorrectForPhylogeny.bf')

        if not exists(batchfile):
            raise ValueError('pass a valid (and existent) batchfile to run')

        fd, self.__inputfile = mkstemp(); close(fd)

        self.__batchfile = batchfile
        self.__commands = commands
        self.__hyphy = HyPhy()

    def __del__(self):
        for file in (self.__inputfile,):
            if file and exists(file):
                remove(file)

    def run(self, seqrecords):
        with open(self.__inputfile, 'w') as fh:
            SeqIO.write(seqrecords, fh, 'fasta') 

        self.__hyphy.ExecuteBF('ExecuteAFile("%s", { "0": "%s" })' % (self.__batchfile, self.__inputfile))

        order = self.__hyphy.AskFor('order').castToString().split(',')
        mat = self.__hyphy.AskFor('data').castToMatrix()

        ret = zeros((mat.mRows, mat.mCols), dtype=float)

        for i in xrange(mat.mRows):
            for j in xrange(mat.mCols):
                ret[i, j] = mat.MatrixCell(i, j)

        return order, ret
