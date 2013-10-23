#!/usr/bin/env python

from __future__ import division, print_function

from argparse import ArgumentParser
from sys import argv as sys_argv, exit as sys_exit

from idepi.encoder import AminoEncoder
from idepi.argument import PathType
from idepi.feature_extraction import MSAVectorizer
from idepi.phylogeny import (
    Phylo,
    PhyloGzFile
)
from idepi.util import (
    is_refseq,
    load_stockholm,
    reference_index
)


def main(args=None):
    if args is None:
        args = sys_argv[1:]

    parser = ArgumentParser(description='Generate a phylogeny from an alignment.')
    parser.add_argument('ALIGNMENT', type=PathType)
    parser.add_argument('OUTPUT', type=PathType)

    ns = parser.parse_args(args)

    msa = load_stockholm(ns.ALIGNMENT)

    try:
        refidx = reference_index(msa, is_refseq)
    except IndexError:
        raise RuntimeError('No reference sequence found!')

    labels = MSAVectorizer(AminoEncoder).fit(msa).get_feature_names()

    seqrecords = [r for i, r in enumerate(msa) if not i == refidx]

    tree, alignment = Phylo()(seqrecords)

    PhyloGzFile.write(ns.OUTPUT, tree, alignment, labels)

    return 0

if __name__ == '__main__':
    sys_exit(main())
