
from contextlib import closing
from multiprocessing import cpu_count
from os import close, remove
from os.path import abspath, exists, join, split
from sys import stderr
from tempfile import mkstemp

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from Bio import SeqIO

from hypy import HyphyInterface


__all__ = ['Ancestors']


class Ancestors(HyphyInterface):

    def __init__(self, batchfile=None, num_cpus=None):
        if batchfile is None:
            batchfile = join(
                split(abspath(__file__))[0],
                'hyphy', 'GetAncestors.bf'
            )

        if num_cpus is None:
            num_cpus = cpu_count()

        super(Ancestors, self).__init__(batchfile, num_cpus)

        self.__inputfile = None

    def __del__(self):
        if self.__inputfile is not None and exists(self.__inputfile):
            remove(self.__inputfile)

    def ancestors(self, seqs, quiet=True):

        if self.__inputfile is None or not exists(self.__inputfile):
            fd, self.__inputfile = mkstemp(); close(fd)

        with open(self.__inputfile, 'w') as fh:
            SeqIO.write(seqs, fh, 'fasta')

        self.queuevar('_inputFile', self.__inputfile)
        self.runqueue()

        if not quiet:
            if self.stdout != '':
                print >> stderr, self.stdout
            if self.warnings != '':
                print >> stderr, self.warnings

        if self.stderr != '':
            raise RuntimeError(self.stderr)

        tree = self.getvar('tree', HyphyInterface.STRING)
        with closing(StringIO(self.getvar('ancestors', HyphyInterface.STRING))) as fh:
            fh.seek(0)
            ancestors = [r for r in SeqIO.parse(fh, 'fasta')]

        if tree[-1] != ';':
            tree += ';'

        return tree, ancestors

    def __call__(self, seqs, quiet=True):
        return Ancestors.ancestors(self, seqs, quiet)
