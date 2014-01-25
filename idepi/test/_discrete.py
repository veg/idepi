
import sys

from os import close, remove
from tempfile import mkstemp

import numpy as np

from BioExt.io import LazyAlignIO as AlignIO

from sklearn.svm import SVC

from sklmrmr import MRMR

from idepi.encoder import AminoEncoder, StanfelEncoder
from idepi.feature_extraction import MSAVectorizer
from idepi.labeledmsa import LabeledMSA
from idepi.labeler import Labeler
from idepi.util import (
    reference_index,
    is_refseq,
    seqrecord_get_values,
    set_util_params
)

from ._common import (
    TEST_AMINO_STO,
    TEST_AMINO_NAMES,
    TEST_STANFEL_NAMES,
    TEST_Y,
    TEST_AMINO_X,
    TEST_STANFEL_X
)


__all__ = ['test_discrete']


def test_discrete(ARGS):
    # set these to this so we don't exclude anything (just testing file generation and parsing)
    ARGS.NUM_FEATURES = 15 # should be enough, the number is known to be 13
    ARGS.MRMR_METHOD = 'MID'
    ARGS.MAX_CONSERVATION = 1.0
    ARGS.MAX_GAP_RATIO    = 1.0
    ARGS.MIN_CONSERVATION = 1.0
    ARGS.CUTOFF = 20.

    # if we don't do this, DOOMBUNNIES
    set_util_params(ARGS.REFSEQ_IDS)

    fd, sto_filename = mkstemp(); close(fd)

    try:
        fh = open(sto_filename, 'w')
        print(TEST_AMINO_STO, file=fh)
        fh.close()

        alignment = AlignIO.read(sto_filename, 'stockholm')

        for ARGS.ENCODER in (AminoEncoder, StanfelEncoder):

            if ARGS.ENCODER == StanfelEncoder:
                TEST_NAMES = TEST_STANFEL_NAMES
                TEST_X = TEST_STANFEL_X
            else:
                TEST_NAMES = TEST_AMINO_NAMES
                TEST_X = TEST_AMINO_X

            # test mRMR and LSVM file generation
            ylabeler = Labeler(
                seqrecord_get_values,
                lambda row: is_refseq(row) or False, # TODO: again filtration function
                lambda x: x > ARGS.CUTOFF,
                False
            )
            alignment, y, ic50 = ylabeler(alignment)

            refidx = reference_index(alignment, is_refseq)
            alignment = LabeledMSA.from_msa_with_ref(alignment, refidx)
            extractor = MSAVectorizer(ARGS.ENCODER)
            x = extractor.fit_transform(alignment)
            colnames = extractor.get_feature_names()

            # test the feature names portion
            try:
                assert(len(colnames) == len(TEST_NAMES))
            except AssertionError:
                raise AssertionError('gen:   %s\ntruth: %s' % (colnames, TEST_NAMES))

            for name in TEST_NAMES:
                try:
                    assert(name in colnames)
                except AssertionError:
                    raise AssertionError('ERROR: \'%s\' not found in %s' % (name, ', '.join(colnames)))

            assert(np.all(TEST_X == x))

            assert(np.all(TEST_Y == y))

            # generate and test the mRMR portion
            mrmr = MRMR(
                estimator=SVC(kernel='linear'),
                n_features_to_select=ARGS.NUM_FEATURES,
                method=ARGS.MRMR_METHOD,
                normalize=ARGS.MRMR_NORMALIZE,
                similar=ARGS.SIMILAR
                )

            mrmr.fit(x, y)

    finally:
        remove(sto_filename)

    print('ALL TESTS PASS', file=sys.stderr)
