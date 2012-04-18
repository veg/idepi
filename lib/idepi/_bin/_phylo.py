#!/usr/bin/env python

from __future__ import division, print_function

from gzip import GzipFile
from os.path import basename, exists
from re import compile as re_compile
from six import StringIO
from sys import argv as sys_argv, exit as sys_exit, stderr, stdout


from idepi import Alphabet, Phylo, PhyloGzFile, column_labels, crude_sto_read, is_HXB2


def main(argv=sys_argv):
    name = basename(argv.pop(0))

    try:
        assert(len(argv) == 2 and exists(argv[0]))
    except:
        print('usage: %s ALIGNMENT OUTPUT' % name, file=stderr)
        sys_exit(-1)

    msa, refseq_offs = crude_sto_read(argv[0], is_HXB2, description=True)

    try:
        refseq = [r for r in msa if is_HXB2(r)][0]
    except IndexError:
        raise RuntimeError('No reference sequence found!')

    seqrecords = [r for r in msa if not is_HXB2(r)]

    tree, alignment = Phylo()(seqrecords)

    PhyloGzFile.write(argv[1], tree, alignment, column_labels(refseq, refseq_offs))

    return 0

if __name__ == '__main__':
    sys_exit(main())
