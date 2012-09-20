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

from math import copysign
from os import close, remove
from os.path import basename, exists, splitext
from re import compile as re_compile
from tempfile import mkstemp

import numpy as np

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq

from idepi import (
    __version__
)
from idepi.alphabet import Alphabet
from idepi.argument import (
    PathType,
    init_args,
    hmmer_args,
    featsel_args,
    mrmr_args,
    optstat_args,
    encoding_args,
    filter_args,
    svm_args,
    cv_args,
    parse_args,
    finalize_args,
    abtypefactory
)
from idepi.databuilder import DataBuilder
from idepi.filter import naivefilter
from idepi.hmmer import Hmmer
from idepi.labeler import Labeler
from idepi.phylogeny import Phylo, PhyloGzFile
from idepi.results import IdepiResults
from idepi.svm import LinearSvm
from idepi.test import test_discrete
from idepi.util import (
    alignment_identify_refidx,
    generate_alignment,
    is_refseq,
    C_range,
    set_util_params,
    seqrecord_get_ic50s,
    seqrecord_set_ic50
)

from mrmr import DiscreteMrmr, PhyloMrmr

from pyxval import CrossValidator, DiscretePerfStats, SelectingGridSearcher


def main(args=None):
    np.seterr(all='raise')

    if args is None:
        args = sys.argv[1:]

    parser, ns, args = init_args(description='Predict IC50 for unlabeled sequences.', args=args)

    parser = hmmer_args(parser)
    parser = featsel_args(parser)
    parser = mrmr_args(parser)
    parser = optstat_args(parser)
    parser = encoding_args(parser)
    parser = filter_args(parser)
    parser = svm_args(parser)
    parser = cv_args(parser)

    parser.add_argument('--tree', dest='TREE')
    parser.add_argument('ANTIBODY', type=abtypefactory(ns.DATA))
    parser.add_argument('SEQUENCES', type=PathType)

    ARGS = parse_args(parser, args, namespace=ns)

    # do some argument parsing
    if ARGS.TEST:
        test_discrete(ARGS)
        finalize_args(ARGS)
        return {}

    similar = True
    if ARGS.MRMR_METHOD == DiscreteMrmr.MAXREL:
        similar = False

    # set the util params
    set_util_params(ARGS.REFSEQ_IDS, ARGS.IC50)

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=ARGS.ALPHABET)

    ab_basename = ''.join((
        ARGS.ANTIBODY,
        '_dna' if ARGS.DNA else '_amino',
        '_clonal' if ARGS.CLONAL else ''
    ))
    alignment_basename = '_'.join((
        ab_basename,
        ARGS.DATA.basename_root,
        __version__
    ))
    fasta_basename = '_'.join((
        ab_basename,
        ARGS.DATA.basename_root,
        splitext(basename(ARGS.SEQUENCES))[0],
        __version__
    ))

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal = ARGS.DATA.seqrecords(ARGS.ANTIBODY, ARGS.CLONAL, ARGS.DNA)

    # if clonal isn't supported, fallback to default
    if clonal != ARGS.CLONAL:
        ab_basename = ''.join(ab_basename.rsplit('_clonal', 1))
        alignment_basename = ''.join(alignment_basename.rsplit('_clonal', 1))

    alignment = generate_alignment(seqrecords, alignment_basename, is_refseq, ARGS)
    refidx = alignment_identify_refidx(alignment, is_refseq)

    seqfilefmt = 'stockholm' if splitext(ARGS.SEQUENCES)[1].find('sto') == 1 else 'fasta'

    # create a temporary file wherein space characters have been removed
    try:
        fd, tmpseq = mkstemp(); close(fd)
        with open(ARGS.SEQUENCES) as oldfh:
            with open(tmpseq, 'w') as tmpfh:
                re_spaces = re_compile(r'[-._]+')
                def spaceless_seqrecord(record):
                    # you must clear the letter annotations before altering seq or BioPython barfs
                    record.letter_annotations.clear()
                    record.seq = Seq(re_spaces.sub('', str(record.seq)), record.seq.alphabet)
                    return record
                SeqIO.write([spaceless_seqrecord(record) for record in SeqIO.parse(oldfh, seqfilefmt)], tmpfh, 'fasta')

        fasta_stofile = fasta_basename + '.sto'
        if not exists(fasta_stofile):
            hmmer = Hmmer(ARGS.HMMER_ALIGN_BIN, ARGS.HMMER_BUILD_BIN)
            hmmer.align(
                alignment_basename + '.hmm',
                tmpseq,
                output=fasta_stofile,
                alphabet=Hmmer.DNA if ARGS.DNA else Hmmer.AMINO,
                outformat=Hmmer.PFAM
            )
    finally:
        if exists(tmpseq):
            remove(tmpseq)

    with open(fasta_stofile) as fh:
        fasta_aln = AlignIO.read(fh, 'stockholm')

    try:
        assert(alignment.get_alignment_length() == fasta_aln.get_alignment_length())
    except AssertionError:
        print(alignment.get_alignment_length(), fasta_aln.get_alignment_length())
        raise

    filter = naivefilter(
        ARGS.MAX_CONSERVATION,
        ARGS.MIN_CONSERVATION,
        ARGS.MAX_GAP_RATIO
    )
    builder = DataBuilder(
        alignment,
        alph,
        refidx,
        filter
    )
    xt = builder(alignment, refidx)
    colnames = builder.labels

    if ARGS.TREE is not None:
        tree, fasta_aln = Phylo()([r for r in fasta_aln if not is_refseq(r)])
    xp = builder(fasta_aln)

    # compute features
    ylabeler = Labeler(
        seqrecord_get_ic50s,
        lambda row: is_refseq(row) or False, # TODO: again filtration function
        lambda x: x > ARGS.IC50GT,
        ARGS.AUTOBALANCE
    )
    yt, ic50 = ylabeler(alignment)
    assert(
        (ic50 is None and not ARGS.AUTOBALANCE) or
        (ic50 is not None and ARGS.AUTOBALANCE)
    )
    if ARGS.AUTOBALANCE:
        ARGS.IC50 = ic50

    if ARGS.MRMR_NORMALIZE:
        DiscreteMrmr._NORMALIZED = True

    sgs = SelectingGridSearcher(
        classifier_cls=LinearSvm,
        selector_cls=PhyloMrmr if ARGS.PHYLOFILTER else DiscreteMrmr,
        validator_cls=CrossValidator,
        gridsearch_kwargs={ 'C': C_range(*ARGS.LOG2C) },
        classifier_kwargs={},
        selector_kwargs={
            'num_features': ARGS.NUM_FEATURES,
            'method': ARGS.MRMR_METHOD
        },
        validator_kwargs={
            'folds': ARGS.CV_FOLDS,
            'scorer_cls': DiscretePerfStats,
            'scorer_kwargs': { 'optstat': ARGS.OPTSTAT }
        },
        weights_func='weights' # we MUST specify this or it will be set to lambda: None
    )

    sgs.learn(xt, yt)
    yp = sgs.predict(xp).astype(int)

    sgs_weights = sgs.weights()
    sgs_features = sgs.features()

    weights = [{
        'position': colnames[featidx],
        'value': int(copysign(1, sgs_weights[idx])) if idx < len(sgs_weights) else 0
    } for idx, featidx in enumerate(sgs_features)]

    ret = IdepiResults(similar)
    ret.metadata(ARGS, len(alignment)-1, np.mean(yt), ARGS.ANTIBODY, forward_select=None)
    ret.weights(weights)
    ret.predictions([row.id for i, row in enumerate(fasta_aln) if i != refidx], yp)

    print(ret.dumps(), file=ARGS.OUTPUT)

    if ARGS.TREE is not None:
        tree_seqrecords = [seqrecord_set_ic50(r, yp[i]) for i, r in enumerate(fasta_aln)]
        PhyloGzFile.write(ARGS.TREE, tree, tree_seqrecords, colnames, ret['metadata'])

    finalize_args(ARGS)

    return ret


if __name__ == '__main__':
    sys.exit(main())
