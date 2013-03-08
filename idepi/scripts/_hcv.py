
import sys

from argparse import ArgumentParser, FileType
from math import copysign
from os import getenv
from random import seed
from re import compile as re_compile

from BioExt.io import LazyAlignIO as AlignIO

from idepi.alphabet import Alphabet
from idepi.argument import (
    cv_args,
    filter_args,
    optstat_args,
    svm_args
    )
from idepi.databuilder import (
    DataBuilder,
    DataBuilderPairwise,
    DataBuilderRegex,
    DataBuilderRegexPairwise,
    DataReducer
    )
from idepi.filter import naivefilter
from idepi.globber import RegexGlobber
from idepi.labeler import Labeler
from idepi.results import Results
from idepi.scorer import Scorer
from idepi.util import C_range

from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFE
from sklearn.grid_search import GridSearchCV
from sklearn.preprocessing import scale
from sklearn.svm import SVC


def is_refseq(seq):
    return 'gi' in seq.id and '329739' in seq.id

def label(seq):
    x = seq.id.split('_')[0].rstrip('1A').rstrip('1B').upper()
    if x == 'HCC':
        return 1
    elif x == 'CHC':
        return -1
    elif x == 'CIRR':
        return -1
    else:
        return None

def main(args=None):

    if args is None:
        args = sys.argv[1:]

    parser = ArgumentParser(description='differentiate HCC, CHC, and cirrhosis')

    parser = cv_args(parser)
    parser = filter_args(parser)
    parser = optstat_args(parser)
    parser = svm_args(parser)

    parser.add_argument('--numfeats', type=int, dest='NUM_FEATURES', default=10)
    parser.add_argument('--step', type=int, dest='STEP', default=1)

    parser.add_argument('--seed', type=int, dest='SEED', default=42)

    # the default of 10% is too low for this analysis
    parser.set_defaults(MAX_GAP_RATIO=1.0)

    parser.add_argument('ALIGNMENT', type=FileType('r'))

    ARGS = parser.parse_args(args)

    seed(ARGS.SEED)

    alphabet = Alphabet(Alphabet.DNA)

    alignment = AlignIO.read(ARGS.ALIGNMENT, 'fasta')

    # do not break, as we need to completely iterate through alignment
    # so as to hit StopIteration and reset the __iter__ generator
    refidx = -1
    for i, seq in enumerate(alignment):
        if is_refseq(seq):
            refidx = i

    if refidx < 0:
        raise RuntimeError('no reference found in alignment')

    filter_ = naivefilter(
        ARGS.MAX_CONSERVATION,
        ARGS.MIN_CONSERVATION,
        ARGS.MAX_GAP_RATIO
        )

    globber = RegexGlobber(
        [re_compile(r'^[^_]+_([^_]+)_HAP[0-9]+_(?P<weight>[0-9]\.[0-9]+(?:e-[0-9]+)?)')],
        alignment,
        refidx
        )

    ylabeler = Labeler(
        label,
        is_refseq,
        lambda x: x
        )
    alignment, y, _ = ylabeler(alignment, globber=globber)

    # do not break, as we need to completely iterate through alignment
    # so as to hit StopIteration and reset the __iter__ generator
    refidx = -1
    for i, seq in enumerate(alignment):
        if is_refseq(seq):
            refidx = i

    if refidx < 0:
        raise RuntimeError('no reference found in alignment')

    builder = DataReducer(
        DataBuilder(
            alignment,
            alphabet,
            refidx,
            filter_
            ),
        DataBuilderPairwise(
            alignment,
            alphabet,
            refidx,
            filter_,
            ARGS.RADIUS
            ) # ,
#         DataBuilderRegex(
#
#             ),
#         DataBuilderRegexPairwise(
#
#             )
        )
    X = builder(
        alignment,
        refidx,
        globber=globber
        # normalize=True
        )

    # scale the data to 0 mean and unit variance
    X = scale(X)

    scorer = Scorer(ARGS.OPTSTAT)
    results = Results(builder.labels, scorer)

    if ARGS.CV_FOLDS > len(y):
        ARGS.CV_FOLDS = len(y)

    for train_idxs, test_idxs in StratifiedKFold(y, ARGS.CV_FOLDS):

        if train_idxs.sum() < 1 or test_idxs.sum() < 1:
            continue

        rfe = RFE(
            estimator=SVC(kernel='linear'),
            n_features_to_select=ARGS.NUM_FEATURES,
            step=ARGS.STEP
            )

        clf = GridSearchCV(
            estimator=rfe,
            param_grid={
                'estimator_params': [dict(C=c) for c in C_range(*ARGS.LOG2C)]
                },
            score_func=scorer,
            n_jobs=int(getenv('NCPU', -1)), # use all but 1 cpu
            pre_dispatch='2 * n_jobs'
            )

        X_train = X[train_idxs]
        y_train = y[train_idxs]

        clf.fit(X_train, y_train)

        X_test = X[test_idxs]
        y_true = y[test_idxs]
        y_pred = clf.predict(X_test)

        coefs = {}
        c = 0
        for i, selected in enumerate(clf.best_estimator_.support_):
            if selected:
                coefs[i] = int(
                    copysign(1, clf.best_estimator_.estimator_.coef_[0, c])
                    )
                c += 1

        results.add(y_true, y_pred, coefs)

    print(results.dumps())

    if ARGS.ALIGNMENT != sys.stdin:
        ARGS.ALIGNMENT.close()

    return 0
