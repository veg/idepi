#!/usr/bin/env python2.7

from gzip import GzipFile
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
        print >> stderr, 'usage: %s ALIGNMENT OUTPUT' % name
        sys_exit(-1)

    with open(argv[0]) as fh:
        seqrecords = [r for r in SeqIO.parse(fh, 'stockholm')]

    refseq = None
    for i, r in enumerate(seqrecords):
        if is_HXB2(r.id):
            refseq = str(seqrecords[i].seq).upper()
            del seqrecords[i]
            break

    if refseq is None:
        raise RuntimeError('no reference sequence was found in the alignment, aborting')

    tree, alignment = Phylo()(seqrecords)

    with GzipFile(argv[1], 'w') as fh:
        print >> fh, '\n'.join(['BEGIN NEWICK', tree, 'END NEWICK', 'BEGIN FASTA'])
        SeqIO.write(alignment, fh, 'fasta')
        print >> fh, '\n'.join(['END FASTA', 'BEGIN REFSEQ', refseq, 'END REFSEQ'])

    return 0

if __name__ == '__main__':
    sys_exit(main())