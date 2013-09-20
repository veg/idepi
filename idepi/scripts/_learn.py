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

from functools import partial
from math import copysign
from os import close, getenv, remove
from os.path import exists, splitext
from re import compile as re_compile, I as re_I
from tempfile import mkstemp

import numpy as np

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

from BioExt.misc import gapless, translate

from idepi import (
    __version__
)
from idepi.alphabet import Alphabet
from idepi.argument import (
    PathType,
    init_args,
    hmmer_args,
    featsel_args,
    feature_args,
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
from idepi.databuilder import (
    DataBuilder,
    DataBuilderPairwise,
    DataBuilderRegex,
    DataBuilderRegexPairwise,
    DataReducer
    )
from idepi.filter import naivefilter
from idepi.hmmer import HMMER
from idepi.labeler import (
    Labeler,
    expression
    )
from idepi.phylogeny import Phylo, PhyloGzFile
from idepi.results import Results
from idepi.scorer import Scorer
from idepi.test import test_discrete
from idepi.util import (
    alignment_identify_refidx,
    generate_alignment,
    is_refseq,
    C_range,
    set_util_params,
    seqrecord_set_values
)

from sklearn.grid_search import GridSearchCV
from sklearn.svm import SVC

from sklmrmr import MRMR


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    np.seterr(all='raise')

    parser, ns, args = init_args(description='Predict IC50 for unlabeled sequences.', args=args)

    parser = hmmer_args(parser)
    parser = featsel_args(parser)
    parser = feature_args(parser)
    parser = mrmr_args(parser)
    parser = optstat_args(parser)
    parser = encoding_args(parser)
    parser = filter_args(parser)
    parser = svm_args(parser)
    parser = cv_args(parser)

    parser.add_argument('--tree', dest='TREE')
    parser.add_argument('ANTIBODY', type=abtypefactory(ns.DATA), nargs='+')
    parser.add_argument('SEQUENCES', type=PathType)

    ARGS = parse_args(parser, args, namespace=ns)

    antibodies = tuple(ARGS.ANTIBODY)

    # convert the hxb2 reference to amino acid, including loop definitions
    if not ARGS.DNA:
        ARGS.REFSEQ = translate(ARGS.REFSEQ)

    # do some argument parsing
    if ARGS.TEST:
        test_discrete(ARGS)
        finalize_args(ARGS)
        return {}

    if ARGS.MRMR_METHOD == 'MAXREL':
        ARGS.SIMILAR = 0.0

    # set the util params
    set_util_params(ARGS.REFSEQ_IDS)

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=ARGS.ALPHABET)

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal, antibodies = ARGS.DATA.seqrecords(antibodies, ARGS.CLONAL, ARGS.DNA)

    ab_basename = ''.join((
        '+'.join(antibodies),
        '_dna' if ARGS.DNA else '_amino',
        '_clonal' if clonal else ''
        ))
    alignment_basename = '_'.join((
        ab_basename,
        ARGS.DATA.basename_root,
        __version__
        ))

    generate_alignment(seqrecords, alignment_basename, is_refseq, ARGS, load=False)

    # remember the number of original sequences
    seqrecords.append(gapless(ARGS.REFSEQ))
    partition = len(seqrecords)

    with open(ARGS.SEQUENCES) as seqfh:
        seqfmt = 'stockholm' if splitext(ARGS.SEQUENCES)[1].find('sto') == 1 else 'fasta'
        seqrecords.extend(gapless(record) for record in SeqIO.parse(seqfh, seqfmt))

     # create a temporary file wherein space characters have been removed
    try:
        fd, tmpseq = mkstemp(); close(fd)
        fd, tmpaln = mkstemp(); close(fd)

        with open(tmpseq, 'w') as tmpfh:
            SeqIO.write((HMMER.valid(record) for record in seqrecords), tmpfh, 'fasta')

        if not exists(alignment_basename + '.hmm'):
            raise RuntimeError('missing HMM profile for alignment')

        hmmer = HMMER(ARGS.HMMER_ALIGN_BIN, ARGS.HMMER_BUILD_BIN)
        hmmer.align(
            alignment_basename + '.hmm',
            tmpseq,
            output=tmpaln,
            alphabet=HMMER.DNA if ARGS.DNA else HMMER.AMINO,
            outformat=HMMER.PFAM
            )

        if not exists(tmpaln):
            raise RuntimeError('unable to align test sequences')

        with open(tmpaln) as tmpfh:
            alignment = AlignIO.read(tmpfh, 'stockholm')
    except ValueError:
        with open(tmpaln) as tmpfh:
            print(tmpfh.read(), file=sys.stderr)
        raise
    finally:
        if exists(tmpseq):
            remove(tmpseq)
        if exists(tmpaln):
            remove(tmpaln)

    train_msa = MultipleSeqAlignment([], alphabet=alignment._alphabet)
    test_msa = MultipleSeqAlignment([], alphabet=alignment._alphabet)
    for i, record in enumerate(alignment):
        if i < partition:
            train_msa.append(record)
        else:
            test_msa.append(record)

    re_pngs = re_compile(r'N[^P][TS][^P]', re_I)

    # this is stupid but works generally speaking
    if ARGS.AUTOBALANCE:
        ARGS.LABEL = ARGS.LABEL.split()[0]

    # compute features
    ylabeler = Labeler(
        partial(expression, label=ARGS.LABEL),
        is_refseq,  # TODO: again, filtration function
        ARGS.AUTOBALANCE
    )
    train_msa, y_train, threshold = ylabeler(train_msa)
    assert (
        (threshold is None and not ARGS.AUTOBALANCE) or
        (threshold is not None and ARGS.AUTOBALANCE)
        )
    if ARGS.AUTOBALANCE:
        ARGS.LABEL = '{0} > {1}'.format(ARGS.LABEL.strip(), threshold)

    filter = naivefilter(
        ARGS.MAX_CONSERVATION,
        ARGS.MIN_CONSERVATION,
        ARGS.MAX_GAP_RATIO
    )

    refidx = alignment_identify_refidx(train_msa, is_refseq)

    builders = [
        DataBuilder(
            alignment,
            alph,
            refidx,
            filter
            )
        ]

    if ARGS.RADIUS:
        builders.append(
            DataBuilderPairwise(
                alignment,
                alph,
                refidx,
                filter,
                ARGS.RADIUS
                )
            )

    if ARGS.PNGS:
        builders.append(
            DataBuilderRegex(
                alignment,
                alph,
                refidx,
                re_pngs,
                4,
                label='PNGS'
                )
            )

    if ARGS.PNGS_PAIRS:
        builders.append(
            DataBuilderRegexPairwise(
                alignment,
                alph,
                refidx,
                re_pngs,
                4,
                label='PNGS'
                )
            )

    builder = DataReducer(*builders)
    X_train = builder(train_msa, refidx)

    if ARGS.TREE is not None:
        tree, test_msa = Phylo()(r for r in test_msa if not is_refseq(r))

    X_pred = builder(test_msa)

    svm = SVC(kernel='linear')

    mrmr = MRMR(
        estimator=svm,
        n_features_to_select=ARGS.NUM_FEATURES,
        method=ARGS.MRMR_METHOD,
        normalize=ARGS.MRMR_NORMALIZE,
        similar=ARGS.SIMILAR
        )

    mrmr.fit(X_train, y_train)

    # train only using the MRMR-selected features
    X_train_ = X_train[:, mrmr.support_]

    scorer = Scorer(ARGS.OPTSTAT)

    clf = GridSearchCV(
        estimator=svm,
        param_grid={'C': list(C_range(*ARGS.LOG2C))},
        scoring=scorer,
        n_jobs=int(getenv('NCPU', -1)),  # use all but 1 cpu
        pre_dispatch='2 * n_jobs',
        cv=ARGS.CV_FOLDS
        )

    clf.fit(X_train_, y_train)

    coef_ = clf.best_estimator_.coef_
    ranking_ = mrmr.ranking_
    support_ = mrmr.support_
    # use only the MRMR-selected features, like above
    X_pred_ = X_pred[:, support_]
    y_pred = clf.predict(X_pred_)

    coefs = {}
    ranks = {}
    col = 0
    for (i, (rank, selected)) in enumerate(zip(ranking_, support_)):
        if selected:
            coefs[i] = int(
                copysign(1, coef_[0, col])
                )
            ranks[i] = int(rank)
            col += 1

    results = Results(builder.labels, scorer, ARGS.SIMILAR)
    results.add(y_train, None, coefs, ranks)
    results.metadata(antibodies, ARGS.LABEL)
    results.predictions([r.id for r in test_msa if not is_refseq(r)], y_pred)

    print(results.dumps(), file=ARGS.OUTPUT)

    if ARGS.TREE is not None:
        tree_seqrecords = [
            seqrecord_set_values(r, 'IC50', [y_pred[i]])
            for i, r
            in enumerate(test_msa)
            ]
        PhyloGzFile.write(ARGS.TREE, tree, tree_seqrecords, builder.labels, results['metadata'])

    finalize_args(ARGS)

    return results


if __name__ == '__main__':
    sys.exit(main())
