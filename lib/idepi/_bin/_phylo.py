#!/usr/bin/env python2.7

from gzip import GzipFile
from io import StringIO
from os.path import basename, exists
from re import compile as re_compile
from sys import argv as sys_argv, exit as sys_exit, stderr, stdout


from Bio import SeqIO

from idepi import Phylo, is_HXB2


def main(argv=sys_argv):
    name = basename(argv.pop(0))

    try:
        assert(len(argv) == 2 and exists(argv[0]))
    except:
        print('usage: %s ALIGNMENT OUTPUT' % name, file=stderr)
        sys_exit(-1)

    with open(argv[0]) as fh:
        seqrecords = [r for r in SeqIO.parse(fh, 'stockholm')]

    refseq = None
    for i, r in enumerate(seqrecords):
        if is_HXB2(r):
            refseq = str(seqrecords[i].seq).upper()
            del seqrecords[i]
            break

    if refseq is None:
        raise RuntimeError('no reference sequence was found in the alignment, aborting')

    tree, alignment = Phylo()(seqrecords)

    with GzipFile(argv[1], 'wb') as fh:
        fh.write(bytes('\n'.join(['BEGIN NEWICK', tree, 'END NEWICK', 'BEGIN FASTA', '']), 'utf-8'))
        with StringIO() as fh2:
            SeqIO.write(alignment, fh2, 'fasta')
            fh.write(bytes(fh2.getvalue().strip(), 'utf-8'))
        fh.write(bytes('\n'.join(['', 'END FASTA', 'BEGIN REFSEQ', refseq, 'END REFSEQ']), 'utf-8'))

    return 0

if __name__ == '__main__':
    sys_exit(main())
