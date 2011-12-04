
from contextlib import closing
from multiprocessing import cpu_count
from os import close, remove
from os.path import abspath, exists, join, split
from re import compile as re_compile, I as re_I
from sys import stderr
from tempfile import mkstemp

try:
    from io import StringIO
except ImportError:
    from io import StringIO

from Bio import SeqIO

from hypy import HyphyInterface


__all__ = ['Phylo']


class Phylo(HyphyInterface):

    def __init__(self, batchfile=None, num_cpus=None):
        if batchfile is None:
            batchfile = join(
                split(abspath(__file__))[0],
                'hyphy', 'GetPhylo.bf'
            )

        if num_cpus is None:
            num_cpus = cpu_count()

        super(Phylo, self).__init__(batchfile, num_cpus)

        self.__inputfile = None

    def __del__(self):
        if self.__inputfile is not None and exists(self.__inputfile):
            remove(self.__inputfile)

    def getphylo(self, seqs, quiet=True):

        if self.__inputfile is None or not exists(self.__inputfile):
            fd, self.__inputfile = mkstemp(); close(fd)

        with open(self.__inputfile, 'w') as fh:
            SeqIO.write(seqs, fh, 'fasta')

        ids = {}
        mangle = re_compile(r'[^a-zA-Z0-9]+', re_I)
        for r in seqs:
            newid = mangle.sub('_', r.id).rstrip('_')
            ids[newid] = r.id

        self.queuevar('_inputFile', self.__inputfile)
        self.runqueue()

        if not quiet:
            if self.stdout != '':
                print(self.stdout, file=stderr)
            if self.warnings != '':
                print(self.warnings, file=stderr)

        if self.stderr != '':
            raise RuntimeError(self.stderr)

        tree = self.getvar('tree', HyphyInterface.STRING)
        with closing(StringIO(self.getvar('ancestors', HyphyInterface.STRING))) as fh:
            fh.seek(0)
            ancestors = [r for r in SeqIO.parse(fh, 'fasta')]

        for r in ancestors:
            key = r.id.rstrip('_unknown_description_').rstrip('_')
            if key in ids:
                newid = ids[key]
                tree = tree.replace(r.id, newid)
                r.id = newid

        if tree[-1] != ';':
            tree += ';'

        return tree, ancestors

    def __call__(self, seqs, quiet=True):
        return Phylo.getphylo(self, seqs, quiet)
