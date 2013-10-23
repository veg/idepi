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
from gzip import open as gzip_open
from os import getenv
from pickle import dump as pickle_dump
from re import compile as re_compile, I as re_I

import numpy as np

from idepi import (
    __version__
)
from idepi.argument import (
    init_args,
    hmmer_args,
    featsel_args,
    feature_args,
    mrmr_args,
    rfe_args,
    optstat_args,
    filter_args,
    svm_args,
    cv_args,
    parse_args,
    finalize_args,
    abtypefactory
)
from idepi.encoder import DNAEncoder
from idepi.feature_extraction import (
    FeatureUnion,
    MSAVectorizer,
    MSAVectorizerPairwise,
    MSAVectorizerRegex,
    MSAVectorizerRegexPairwise
    )
from idepi.filters import naive_filter
from idepi.labeler import (
    Labeler,
    expression
    )
from idepi.results import Results
from idepi.scorer import Scorer
from idepi.test import test_discrete
from idepi.util import (
    coefs_ranks,
    generate_alignment,
    is_refseq,
    C_range,
    set_util_params
)

from sklearn.grid_search import GridSearchCV
from sklearn.svm import SVC

from sklmrmr import MRMR


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    np.seterr(all='raise')

    parser, ns, args = init_args(description='learn model for labeled sequences', args=args)

    parser = hmmer_args(parser)
    parser = featsel_args(parser)
    parser = feature_args(parser)
    parser = mrmr_args(parser)
    parser = rfe_args(parser)
    parser = optstat_args(parser)
    parser = filter_args(parser)
    parser = svm_args(parser)
    parser = cv_args(parser)

    parser.add_argument('--tree', dest='TREE')
    parser.add_argument('ANTIBODY', type=abtypefactory(ns.DATA), nargs='+')
    parser.add_argument('MODEL', type=str)

    ARGS = parse_args(parser, args, namespace=ns)

    antibodies = tuple(ARGS.ANTIBODY)

    # do some argument parsing
    if ARGS.TEST:
        test_discrete(ARGS)
        finalize_args(ARGS)
        return {}

    if ARGS.MRMR_METHOD == 'MAXREL':
        ARGS.SIMILAR = 0.0

    # set the util params
    set_util_params(ARGS.REFSEQ.id)

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal, antibodies = ARGS.DATA.seqrecords(antibodies, ARGS.CLONAL)

    ab_basename = ''.join((
        '+'.join(antibodies),
        '_dna' if ARGS.ENCODER == DNAEncoder else '_amino',
        '_clonal' if clonal else ''
        ))
    alignment_basename = '_'.join((
        ab_basename,
        ARGS.DATA.basename_root,
        __version__
        ))
    sto_filename = alignment_basename + '.sto'

    alignment, hmm = generate_alignment(seqrecords, sto_filename, is_refseq, ARGS)

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
    alignment, y, threshold = ylabeler(alignment)

    assert (
        (threshold is None and not ARGS.AUTOBALANCE) or
        (threshold is not None and ARGS.AUTOBALANCE)
        )
    if ARGS.AUTOBALANCE:
        ARGS.LABEL = '{0} > {1}'.format(ARGS.LABEL.strip(), threshold)

    filter = naive_filter(
        max_conservation=ARGS.MAX_CONSERVATION,
        min_conservation=ARGS.MIN_CONSERVATION,
        max_gap_ratio=ARGS.MAX_GAP_RATIO
        )

    extractors = [('site', MSAVectorizer(ARGS.ENCODER, filter))]

    if ARGS.RADIUS:
        extractors.append(('site_pairs', MSAVectorizerPairwise(ARGS.ENCODER, filter, ARGS.RADIUS)))

    if ARGS.PNGS:
        extractors.append(('pngs', MSAVectorizerRegex(re_pngs, 4, name='PNGS')))

    if ARGS.PNGS_PAIRS:
        extractors.append(
            ('pngs_pairs', MSAVectorizerRegexPairwise(re_pngs, 4, name='PNGS'))
            )

    extractor = FeatureUnion(extractors, n_jobs=1)  # n_jobs must be one for now
    X = extractor.fit_transform(alignment)

    svm = SVC(kernel='linear', class_weight='auto')

    mrmr = MRMR(
        estimator=svm,
        n_features_to_select=ARGS.NUM_FEATURES,
        method=ARGS.MRMR_METHOD,
        normalize=ARGS.MRMR_NORMALIZE,
        similar=ARGS.SIMILAR
        )

    mrmr.fit(X, y)

    # train only using the MRMR-selected features
    X_ = X[:, mrmr.support_]

    scorer = Scorer(ARGS.OPTSTAT)

    clf = GridSearchCV(
        estimator=svm,
        param_grid=dict(C=list(C_range(*ARGS.LOG2C))),
        scoring=scorer,
        n_jobs=int(getenv('NCPU', -1)),  # use all but 1 cpu
        pre_dispatch='2 * n_jobs',
        cv=ARGS.CV_FOLDS
        )

    clf.fit(X_, y)

    with gzip_open(ARGS.MODEL, 'wb') as fh:
        pickle_dump((2, ARGS.ENCODER, ARGS.LABEL, hmm, extractor, mrmr, clf), fh)

    coefs, ranks = coefs_ranks(mrmr.ranking_, mrmr.support_, clf.best_estimator_.coef_)
    results = Results(extractor.get_feature_names(), scorer, ARGS.SIMILAR)

    results.add(y, clf.predict(X_), coefs, ranks)
    results.metadata(antibodies, ARGS.LABEL)

    print(results.dumps(), file=ARGS.OUTPUT)

    finalize_args(ARGS)

    return ARGS.MODEL


if __name__ == '__main__':
    sys.exit(main())
