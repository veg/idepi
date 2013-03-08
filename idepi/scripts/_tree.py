#!/usr/bin/env python

from __future__ import division, print_function

from argparse import ArgumentParser, ArgumentTypeError
from re import compile as re_compile, sub as re_sub
from sys import argv as sys_argv, exit as sys_exit

from numpy import mean

from idepi.argument import PathType
from idepi.phylogeny import PhyloGzFile
from idepi.util import seqrecord_get_values


NUMERIC = re_compile(r'[^0-9]+')


def feattype(string):
    try:
        return set(int(NUMERIC.sub('', f)) for f in string.split(','))
    except ValueError:
        raise ArgumentTypeError('FEATURES must be a comma-delimited list of ints')


def main(args=None):
    if args is None:
        args = sys_argv[1:]

    parser = ArgumentParser(description='')
    parser.add_argument('TREE', type=PathType)
    parser.add_argument('FEATURES', type=feattype)
    ns = parser.parse_args(args)

    tree, alignment, colnames, _ = PhyloGzFile.read(ns.TREE)

    icolnames = [(idx, colname) for idx, colname in enumerate(colnames) if int(NUMERIC.sub('', colname)) in ns.FEATURES]

    for r in alignment:
        # labels has length of icolnames plus the ic50
        labels = [None] * (len(icolnames) + 1)
        i = 1
        for idx, colname in icolnames:
            if len(colnames) > 1:
                labels[i] = colname + r.seq[idx]
            else:
                labels[i] = r.seq[idx]
            i += 1
        try:
            labels[0] = '%.3g' % mean(seqrecord_get_values(r))
        except ValueError:
            if not (len(r.id) > 4 and r.id[:4].lower() == 'node'):
                print(r)
            labels.pop(0)
        # include the ':' here to make sure we grab the end of the label
        tree = re_sub(r'([,()])' + r.id + r'(?:_[0-9]+)?:', r'\g<1>' + '_'.join(labels) + ':', tree)

    print(tree)

    return 0

if __name__ == '__main__':
    sys_exit(main())
