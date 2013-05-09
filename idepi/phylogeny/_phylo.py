
from __future__ import division, print_function

from contextlib import closing
from multiprocessing import cpu_count
from os import close, remove
from os.path import abspath, exists, join, split
from re import compile as re_compile, I as re_I
from six import StringIO
from sys import stderr
from tempfile import mkstemp

from Bio import AlignIO, SeqIO

from hppy import HyphyInterface


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

        seqs = list(seqs)

        if self.__inputfile is None or not exists(self.__inputfile):
            fd, self.__inputfile = mkstemp(); close(fd)

        with open(self.__inputfile, 'w') as fh:
            SeqIO.write(seqs, fh, 'fasta')

        newick_mangle = re_compile(r'[()]')

        id_descs = {}
        mangle = re_compile(r'[^a-zA-Z0-9]+', re_I)
        for r in seqs:
            newid = mangle.sub('_', '_'.join((r.id, r.description))).rstrip('_')
            id_descs[newid] = (newick_mangle.sub('_', r.id).strip('_'), r.description)

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
            ancestors = AlignIO.read(fh, 'fasta')

        hyphymangling = re_compile(r'_[0-9]+$')

        for r in ancestors:
            key = r.id.rstrip('_unknown_description_').rstrip('_')
            if key not in id_descs:
                key = hyphymangling.sub('', key)
                if key not in id_descs:
                    continue
            # if the key exists, replace
            oldid, olddesc = id_descs[key]
            tree = tree.replace(r.id, oldid)
            r.id = oldid
            r.description = olddesc

        if tree[-1] != ';':
            tree += ';'

        return tree, ancestors

    def __call__(self, seqs, quiet=True):
        return Phylo.getphylo(self, seqs, quiet)
