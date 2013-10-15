#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# discrete.py :: computes a predictive model for an nAb in IDEPI's database,
# and provides cross-validation performance statistics and HXB2-relative
# feature information for this model.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import division, print_function

import sys

from gzip import open as gzip_open
from os import remove
from os.path import exists
from pickle import load as pickle_load

import numpy as np

from Bio import SeqIO

from BioExt.misc import gapless, translate

from idepi.argument import (
    PathType,
    init_args,
    hmmer_args,
    parse_args,
    finalize_args
)
from idepi.util import (
    generate_alignment_,
    load_stockholm,
    seqfile_format
    )


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    np.seterr(all='raise')

    parser, ns, args = init_args(description='Predict IC50 for unlabeled sequences.', args=args)

    parser = hmmer_args(parser)

    parser.add_argument('MODEL', type=PathType)
    parser.add_argument('SEQUENCES', type=PathType)

    ARGS = parse_args(parser, args, namespace=ns)

    with gzip_open(ARGS.MODEL, 'rb') as fh:
        ARGS.DNA, ARGS.LABEL, builder, mrmr, clf = pickle_load(fh)

    # create a temporary file wherein space characters have been removed
    with open(ARGS.SEQUENCES) as seq_fh:
        seq_fmt = seqfile_format(ARGS.SEQUENCES)

        def seqrecords():
            for record in SeqIO.parse(seq_fh, seq_fmt):
                r = gapless(record)
                if not ARGS.DNA:
                    r = translate(r)
                yield r

        try:
            tmpaln = generate_alignment_(seqrecords(), ARGS)
            alignment = load_stockholm(tmpaln, trim=True)
        finally:
            if exists(tmpaln):
                remove(tmpaln)

    X = builder(alignment)[:, mrmr.support_]
    y = clf.predict(X)

    idlen = max(len(r.id) for r in alignment) + 3

    print('{\n  "predictions": [', file=ARGS.OUTPUT)
    for i, (r, p) in enumerate(zip(alignment, y)):
        if i > 0:
            print(',')
        print(
            '    {{{{ "id": {{0:<{0:d}s}} "value": {{1: d}} }}}}'.format(idlen).format(
                '"{0:s}",'.format(r.id), p),
            file=ARGS.OUTPUT, end='')
    print('\n  ]\n}', file=ARGS.OUTPUT)

    finalize_args(ARGS)

    return 0


if __name__ == '__main__':
    sys.exit(main())
