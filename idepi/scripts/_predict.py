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
from os import close, remove
from os.path import exists
from pickle import load as pickle_load
from tempfile import mkstemp

import numpy as np

from Bio import SeqIO

from BioExt.misc import translate

from idepi.argument import (
    PathType,
    init_args,
    hmmer_args,
    parse_args,
    finalize_args
)
from idepi.constants import AminoAlphabet, DNAAlphabet
from idepi.encoder import DNAEncoder
from idepi.util import (
    generate_alignment_,
    load_stockholm,
    seqfile_format
    )
from idepi.verifier import VerifyError, Verifier


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    np.seterr(all='raise')

    parser, ns, args = init_args(description='Predict label for unlabeled sequences', args=args)

    parser = hmmer_args(parser)

    parser.add_argument('MODEL', type=PathType)
    parser.add_argument('SEQUENCES', type=PathType)

    ARGS = parse_args(parser, args, namespace=ns)

    with gzip_open(ARGS.MODEL, 'rb') as fh:
        try:
            model = pickle_load(fh)
            if model[0] != 4:
                raise ImportError('incompatible model version')
            ARGS.ENCODER, ARGS.LABEL, hmm, extractor, clf = model[1:]
        except ImportError:
            msg = 'your model is not of the appropriate version, please re-learn your model'
            raise RuntimeError(msg)

    # create a temporary file wherein space characters have been removed
    with open(ARGS.SEQUENCES) as seq_fh:

        def seqrecords():
            is_dna = ARGS.ENCODER == DNAEncoder
            seq_fmt = seqfile_format(ARGS.SEQUENCES)
            source = Verifier(SeqIO.parse(seq_fh, seq_fmt), DNAAlphabet)
            try:
                for record in source:
                    yield record if is_dna else translate(record)
            except VerifyError:
                if is_dna:
                    msg = (
                        "your model specifies a DNA encoding "
                        "which is incompatible with protein sequences"
                        )
                    raise RuntimeError(msg)
                source.set_alphabet(AminoAlphabet)
                for record in source:
                    yield record

        try:
            fd, tmphmm = mkstemp(); close(fd)
            with open(tmphmm, 'wb') as hmm_fh:
                hmm_fh.write(hmm)
                # explicitly gc hmm
                hmm = None
            tmpaln = generate_alignment_(seqrecords(), tmphmm, ARGS)
            alignment = load_stockholm(tmpaln, trim=True)
        finally:
            if exists(tmphmm):
                remove(tmphmm)
            if exists(tmpaln):
                remove(tmpaln)

    X = extractor.transform(alignment)
    y = clf.predict(X)

    feature_names = extractor.get_feature_names()
    support = clf.named_steps['mrmr'].support_
    labels = ['"{0:s}"'.format(feature_names[i]) for i, s in enumerate(support) if s]
    emptys = [' ' * (len(label) + 2) for label in labels]
    idlen = max(len(r.id) for r in alignment) + 3

    print('{{\n  "label": "{0:s}",\n  "predictions": ['.format(ARGS.LABEL), file=ARGS.OUTPUT)
    for i, r in enumerate(alignment):
        if i > 0:
            print(',')
        features = ['[ ']
        for j, x in enumerate(X[i, support]):
            if x:
                features.append(labels[j])
                features.append(', ')
            else:
                features.append(emptys[j])
        features.append(' ]')
        # replace the last comma with a space
        idx = None
        for k, f in enumerate(features):
            if f == ', ':
                idx = k
        if idx is None:
            features[0] = features[0].rstrip()
            features[-1] = features[-1].lstrip()
        else:
            features[idx] = ''
        features_ = ''.join(features)
        print(
            '    {{{{ "id": {{0:<{0:d}s}} "value": {{1: d}}, "features": {{2:s}} }}}}'.format(
                idlen).format('"{0:s}",'.format(r.id), y[i], features_),
            file=ARGS.OUTPUT, end='')
    print('\n  ]\n}', file=ARGS.OUTPUT)

    finalize_args(ARGS)

    return 0


if __name__ == '__main__':
    sys.exit(main())
