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
from os import getenv
from os.path import join
from re import compile as re_compile, I as re_I

import numpy as np

from idepi import (
    __path__ as idepi_path,
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
    finalize_args,
    parse_args,
    AntibodyTypeFactory
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
    expression,
    skipper
    )
from idepi.logging import init_log
from idepi.results import Results
from idepi.scorer import Scorer
from idepi.test import test_discrete
from idepi.util import (
    coefs_ranks,
    generate_alignment,
    is_refseq,
    set_util_params,
    C_range
    )

from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFE
from sklearn.grid_search import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC

from sklmrmr import MRMR


RAND_HIV_ENV_PEP_STOFILE = join(idepi_path[0], 'data', 'randhivenvpep_final.sto')


def main(args=None):
    init_log()

    if args is None:
        args = sys.argv[1:]

    np.seterr(all='raise')

    # so some option parsing
    parser, ns, args = init_args(description="Predict epitope sites.", args=args)

    parser = hmmer_args(parser)
    parser = featsel_args(parser)
    parser = feature_args(parser)
    parser = mrmr_args(parser)
    parser = rfe_args(parser)
    parser = optstat_args(parser)
    parser = filter_args(parser)
    parser = svm_args(parser)
    parser = cv_args(parser)

    parser.add_argument('ANTIBODY', type=AntibodyTypeFactory(ns.DATA), nargs='+')

    ARGS = parse_args(parser, args, namespace=ns)

    # do some argument parsing
    if ARGS.TEST:
        test_discrete(ARGS)
        finalize_args(ARGS)
        return {}

    # maxrel doesn't support similar
    if ARGS.MRMR_METHOD == 'MAXREL':
        ARGS.SIMILAR = 0.0

    antibodies = tuple(ARGS.ANTIBODY)

    # set the util params
    set_util_params(ARGS.REFSEQ.id)

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal, antibodies = ARGS.DATA.seqrecords(antibodies, ARGS.CLONAL)

    # if we're doing LOOCV, make sure we set CV_FOLDS appropriately
    if ARGS.LOOCV:
        ARGS.CV_FOLDS = len(seqrecords)

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

    # don't capture the second variable, let it be gc'd
    alignment = generate_alignment(seqrecords, sto_filename, is_refseq, ARGS)[0]

    re_pngs = re_compile(r'N[^P][TS][^P]', re_I)

    ylabeler = Labeler(
        partial(expression, ARGS.LABEL),
        partial(skipper, is_refseq, ARGS.SUBTYPES)
        )
    alignment, y, threshold = ylabeler(alignment)

    filter = naive_filter(
        max_conservation=ARGS.MAX_CONSERVATION,
        min_conservation=ARGS.MIN_CONSERVATION,
        max_gap_ratio=ARGS.MAX_GAP_RATIO
        )

    extractors = [('site_ident', MSAVectorizer(ARGS.ENCODER, filter))]

    if ARGS.RADIUS:
        extractors.append(('pair_ident', MSAVectorizerPairwise(ARGS.ENCODER, filter, ARGS.RADIUS)))

    if ARGS.PNGS:
        extractors.append(('pngs', MSAVectorizerRegex(re_pngs, 4, name='PNGS')))

    if ARGS.PNGS_PAIRS:
        extractors.append(
            ('pngs_pair', MSAVectorizerRegexPairwise(re_pngs, 4, name='PNGS'))
            )

    extractor = FeatureUnion(extractors, n_jobs=1)  # n_jobs must be 1 for now
    X = extractor.fit_transform(alignment)

    assert y.shape[0] == X.shape[0], \
        "number of classes doesn't match the data: %d vs %d" % (y.shape[0], X.shape[0])

    scorer = Scorer(ARGS.OPTSTAT)

    # do grid-search as part of the svm to avoid
    # performing feature selection on every iteration
    # of the grid search, which naturally takes forever
    svm = GridSearchCV(
        estimator=SVC(kernel='linear', class_weight='auto'),
        param_grid=dict(C=list(C_range(*ARGS.LOG2C))),
        scoring=scorer,
        n_jobs=int(getenv('NCPU', -1)),
        pre_dispatch='3 * n_jobs',
        cv=ARGS.CV_FOLDS - 1
        )

    results = None
    for n_features in ARGS.FEATURE_GRID:
        results_ = Results(extractor.get_feature_names(), scorer, ARGS.SIMILAR)

        for train_idxs, test_idxs in StratifiedKFold(y, ARGS.CV_FOLDS):

            if train_idxs.sum() < 1 or test_idxs.sum() < 1:
                y_true = y[test_idxs]
                results_.add(y_true, y_true, {})
                continue

            X_train = X[train_idxs]
            y_train = y[train_idxs]

            if ARGS.RFE:
                clf = RFE(
                    estimator=svm,
                    n_features_to_select=n_features,
                    step=ARGS.RFE_STEP
                    )
            else:
                mrmr = MRMR(
                    k=n_features,
                    method=ARGS.MRMR_METHOD,
                    normalize=ARGS.MRMR_NORMALIZE,
                    similar=ARGS.SIMILAR
                    )
                clf = Pipeline([('mrmr', mrmr), ('svm', svm)])

            clf.fit(X_train, y_train)

            X_test = X[test_idxs]
            y_true = y[test_idxs]

            if ARGS.RFE:
                selector_ = clf
                svm_ = clf.estimator_.best_estimator_
            else:
                selector_ = clf.named_steps['mrmr']
                svm_ = clf.named_steps['svm'].best_estimator_

            y_pred = clf.predict(X_test)

            coefs, ranks = coefs_ranks(selector_.ranking_, selector_.support_, svm_.coef_)

            results_.add(y_true, y_pred, coefs, ranks)

        if results is None or results_ > results:
            results = results_

    # the alignment reflects the number of sequences either naturally
    results.metadata(antibodies, ARGS.LABEL)

    print(results.dumps(), file=ARGS.OUTPUT)

    finalize_args(ARGS)

    return results


if __name__ == '__main__':
    sys.exit(main())
