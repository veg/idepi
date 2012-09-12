#!/usr/bin/env python

from __future__ import division, print_function

from argparse import ArgumentParser
from sys import argv as sys_argv, exit as sys_exit

from idepi.argument import PathType
from idepi.phylogeny import (
    Phylo,
    PhyloGzFile
)
from idepi.util import (
    alignment_identify_refidx,
    site_labels,
    crude_sto_read,
    is_refseq
)


def main(args=None):
    if args is None:
        args = sys_argv[1:]

    parser = ArgumentParser(description='Generate a phylogeny from an alignment.')
    parser.add_argument('ALIGNMENT', type=PathType)
    parser.add_argument('OUTPUT', type=PathType)

    ns = parser.parse_args(args)

    msa, refseq_offs = crude_sto_read(ns.ALIGNMENT, is_refseq, description=True)

    try:
        refidx = alignment_identify_refidx(msa, is_refseq)
        refseq = msa[refidx]
    except IndexError:
        raise RuntimeError('No reference sequence found!')

    seqrecords = [r for i, r in enumerate(msa) if not i == refidx]

    tree, alignment = Phylo()(seqrecords)

    PhyloGzFile.write(ns.OUTPUT, tree, alignment, site_labels(refseq, refseq_offs))

    return 0

if __name__ == '__main__':
    sys_exit(main())
