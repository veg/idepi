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
from os import getenv
from os.path import join
from re import compile as re_compile, I as re_I
from warnings import warn

import numpy as np

from BioExt.misc import translate

from idepi import (
    __path__ as idepi_path,
    __version__
    )
from idepi.alphabet import Alphabet
from idepi.argument import (
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
    finalize_args,
    parse_args,
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
from idepi.labeler import (
    Labeler,
    expression
    )
from idepi.logging import init_log
from idepi.results import Results
from idepi.scorer import Scorer
from idepi.test import test_discrete
from idepi.util import (
    alignment_identify_refidx,
    generate_alignment,
    is_refseq,
    set_util_params,
    C_range
    )

from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFE
from sklearn.grid_search import GridSearchCV
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
    parser = optstat_args(parser)
    parser = encoding_args(parser)
    parser = filter_args(parser)
    parser = svm_args(parser)
    parser = cv_args(parser)

    parser.add_argument('ANTIBODY', type=abtypefactory(ns.DATA), nargs='+')
    parser.add_argument('--rfe', action='store_true', dest='RFE')
    parser.add_argument('--rfestep', type=int, dest='RFE_STEP')

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

    if len(ARGS.FILTER) != 0:
        if ARGS.NUM_FEATURES > len(ARGS.FILTER):
            ARGS.NUM_FEATURES = len(ARGS.FILTER)
            warn('clamping --numfeats to sizeof(--filter) = %d' % ARGS.NUM_FEATURES)

    # convert the hxb2 reference to amino acid, including loop definitions
    if not ARGS.DNA:
        ARGS.REFSEQ = translate(ARGS.REFSEQ)

    # set the util params
    set_util_params(ARGS.REFSEQ_IDS)

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=ARGS.ALPHABET)

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal, antibodies = ARGS.DATA.seqrecords(antibodies, ARGS.CLONAL, ARGS.DNA)

    # if we're doing LOOCV, make sure we set CV_FOLDS appropriately
    if ARGS.LOOCV:
        ARGS.CV_FOLDS = len(seqrecords)

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

    alignment = generate_alignment(seqrecords, alignment_basename, is_refseq, ARGS)
    filter = naivefilter(
        ARGS.MAX_CONSERVATION,
        ARGS.MIN_CONSERVATION,
        ARGS.MAX_GAP_RATIO
        )

    re_pngs = re_compile(r'N[^P][TS][^P]', re_I)

    # this is stupid but works generally speaking
    if ARGS.AUTOBALANCE:
        ARGS.LABEL = ARGS.LABEL.split()[0]

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

    refidx = alignment_identify_refidx(alignment, is_refseq)

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
    X = builder(alignment, refidx)
    assert y.shape[0] == X.shape[0], \
        "number of classes doesn't match the data: %d vs %d" % (y.shape[0], X.shape[0])

    scorer = Scorer(ARGS.OPTSTAT)

    # compute features
    forward_initval = 1 if ARGS.FORWARD_SELECT else ARGS.NUM_FEATURES
    results = None
    for num_features in range(forward_initval, ARGS.NUM_FEATURES + 1):
        results_ = Results(builder.labels, scorer, ARGS.SIMILAR)

        for train_idxs, test_idxs in StratifiedKFold(y, ARGS.CV_FOLDS):

            if train_idxs.sum() < 1 or test_idxs.sum() < 1:
                y_true = y[test_idxs]
                results_.add(y_true, y_true, {})
                continue

            X_train = X[train_idxs]
            y_train = y[train_idxs]

            svm = SVC(kernel='linear')

            if ARGS.RFE:
                X_train_ = X_train
                estimator = RFE(
                    estimator=svm,
                    n_features_to_select=num_features,
                    step=ARGS.RFE_STEP
                    )
                param_grid = {
                    'estimator_params': [dict(C=c) for c in C_range(*ARGS.LOG2C)]
                    }
            else:
                # do MRMR-fitting separately from GridSearchCV to avoid
                # performing MRMR on every iteration of the grid search,
                # which would naturally take forever
                mrmr = MRMR(
                    estimator=svm,
                    n_features_to_select=num_features,
                    method=ARGS.MRMR_METHOD,
                    normalize=ARGS.MRMR_NORMALIZE,
                    similar=ARGS.SIMILAR
                    )
                mrmr.fit(X_train, y_train)
                # train only using the MRMR-selected features
                X_train_ = X_train[:, mrmr.support_]
                estimator = svm
                param_grid = {
                    'C': list(C_range(*ARGS.LOG2C))
                    }

            clf = GridSearchCV(
                estimator=estimator,
                param_grid=param_grid,
                score_func=scorer,
                n_jobs=int(getenv('NCPU', -1)),  # use all but 1 cpu
                pre_dispatch='2 * n_jobs'
                )

            clf.fit(X_train_, y_train, cv=ARGS.CV_FOLDS - 1)

            X_test = X[test_idxs]
            y_true = y[test_idxs]

            if ARGS.RFE:
                coef_ = clf.best_estimator_.estimator_.coef_
                ranking_ = clf.best_estimator_.ranking_
                support_ = clf.best_estimator_.support_
                X_test_ = X_test
            else:
                coef_ = clf.best_estimator_.coef_
                ranking_ = mrmr.ranking_
                support_ = mrmr.support_
                # use only the MRMR-selected features, like above
                X_test_ = X_test[:, support_]

            y_pred = clf.predict(X_test_)

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

            results_.add(y_true, y_pred, coefs, ranks)

        if results is not None and results_ <= results:
            break

        results = results_

    # the alignment reflects the number of sequences either naturally
    results.metadata(antibodies, ARGS.LABEL)

    print(results.dumps(), file=ARGS.OUTPUT)

    finalize_args(ARGS)

    return results


if __name__ == '__main__':
    sys.exit(main())
