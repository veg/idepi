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

from argparse import ArgumentTypeError
from os.path import join
from re import compile as re_compile, I as re_I
from warnings import warn

import numpy as np

from BioExt import translate

from idepi import (
    __path__ as idepi_path,
    __version__
)
from idepi.alphabet import Alphabet
from idepi.argument import (
    init_args,
    hmmer_args,
    featsel_args,
    mrmr_args,
    optstat_args,
    encoding_args,
    filter_args,
    svm_args,
    cv_args,
    simulation_args,
    finalize_args,
    parse_args,
    abtypefactory
)
from idepi.databuilder import (
    DataBuilder,
    DataBuilderPairwise,
    DataBuilderRegex,
    DataBuilderRegexPairwise,
    DataBuilderRegexTriplewise,
    DataReducer
)
from idepi.filter import naivefilter
from idepi.labeler import Labeler
from idepi.logging import init_log
from idepi.results import IdepiResults
from idepi.simulation import (
    DumbSimulation,
    MarkovSimulation,
    Simulation
)
from idepi.svm import LinearSvm
from idepi.test import test_discrete
from idepi.util import (
    alignment_identify_refidx,
    extract_feature_weights,
    extract_feature_weights_similar,
    generate_alignment,
    is_refseq,
    seqrecord_get_ic50s,
    set_util_params,
    C_range
)

from mrmr import DiscreteMrmr

from pyxval import CrossValidator, DiscretePerfStats, SelectingNestedCrossValidator


RAND_HIV_ENV_PEP_STOFILE = join(idepi_path[0], 'data', 'randhivenvpep_final.sto')


def main(args=None):
    init_log()

    np.seterr(all='raise')

    if args is None:
        args = sys.argv[1:]

    # so some option parsing
    parser, ns, args = init_args(description="Predict epitope sites.", args=args)

    parser = hmmer_args(parser)
    parser = featsel_args(parser)
    parser = mrmr_args(parser)
    parser = optstat_args(parser)
    parser = encoding_args(parser)
    parser = filter_args(parser)
    parser = svm_args(parser)
    parser = cv_args(parser)
    parser = simulation_args(parser)

    parser.add_argument('ANTIBODY', type=abtypefactory(ns.DATA))

    ARGS = parse_args(parser, args, namespace=ns)

    # do some argument parsing
    if ARGS.TEST:
        test_discrete(ARGS)
        finalize_args(ARGS)
        return {}

    similar = True
    if ARGS.MRMR_METHOD == DiscreteMrmr.MAXREL:
        similar = False

    # validate the antibody argument, currently a hack exists to make PG9/PG16 work
    # TODO: Fix pg9/16 hax
    antibody = ARGS.ANTIBODY.strip()

    if ARGS.SIM in (Simulation.EPITOPE, Simulation.SEQUENCE) and ARGS.DNA:
        raise ArgumentTypeError('randseq simulation target not compatible with DNA alphabet')

    if len(ARGS.FILTER) != 0:
        if ARGS.NUM_FEATURES > len(ARGS.FILTER):
            ARGS.NUM_FEATURES = len(ARGS.FILTER)
            warn('clamping --numfeats to sizeof(--filter) = %d' % ARGS.NUM_FEATURES)

    # convert the hxb2 reference to amino acid, including loop definitions
    if not ARGS.DNA:
        ARGS.REFSEQ = translate(ARGS.REFSEQ)

    # set the util params
    set_util_params(ARGS.REFSEQ_IDS, ARGS.IC50)

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=ARGS.ALPHABET)

    sim = None
    if ARGS.SIM is not None:
        if ARGS.SIM == Simulation.DUMB:
            sim = DumbSimulation(ARGS.SIM_RUNS, Simulation.EPITOPE, str(ARGS.REFSEQ.seq))
        elif ARGS.SIM in (Simulation.EPITOPE, Simulation.SEQUENCE):
            sim = MarkovSimulation(ARGS.SIM_RUNS, ARGS.SIM, RAND_HIV_ENV_PEP_STOFILE)
        else:
            raise ValueError("Unknown simulation type '%s'" % ARGS.SIM)

    ab_basename = ''.join((
        antibody,
        '_randseq' if sim is not None and sim.mode == Simulation.SEQUENCE else '',
        '_dna' if ARGS.DNA else '_amino',
        '_clonal' if ARGS.CLONAL else ''
    ))
    alignment_basename = '_'.join((
        ab_basename,
        ARGS.DATA.basename_root,
        __version__
    ))

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal = ARGS.DATA.seqrecords(antibody, ARGS.CLONAL, ARGS.DNA)

    # if we're doing LOOCV, make sure we set CV_FOLDS appropriately
    if ARGS.LOOCV:
        ARGS.CV_FOLDS = len(seqrecords)

    # if clonal isn't supported, fallback to default
    if clonal != ARGS.CLONAL:
        ab_basename = ''.join(ab_basename.rsplit('_clonal', 1))
        alignment_basename = ''.join(alignment_basename.rsplit('_clonal', 1))

    alignment = generate_alignment(seqrecords, alignment_basename, is_refseq, ARGS)
    refidx = alignment_identify_refidx(alignment, is_refseq)
    filter = naivefilter(
        ARGS.MAX_CONSERVATION,
        ARGS.MIN_CONSERVATION,
        ARGS.MAX_GAP_RATIO
    )

    re_pngs = re_compile(r'N[^P][TS][^P]', re_I)

    if sim is None:
        builder = DataReducer(
            DataBuilder(
                alignment,
                alph,
                refidx,
                filter
            ),
            DataBuilderPairwise(
                alignment,
                alph,
                refidx,
                filter,
                ARGS.RADIUS
            ),
            DataBuilderRegex(
                alignment,
                alph,
                refidx,
                re_pngs,
                'PNGS'
            ),
            DataBuilderRegexPairwise(
                alignment,
                alph,
                refidx,
                re_pngs,
                'PNGS'
            ) # ,
#             DataBuilderRegexTriplewise(
#                 alignment,
#                 alph,
#                 refidx,
#                 re_pngs,
#                 'PNGS'
#             )
        )
        x = builder(alignment, refidx)
        colnames = builder.labels
    else:
        if ARGS.SIM_EPI_N is None:
            ARGS.SIM_EPI_N = len(seqrecords)

    # compute features
    forward_initval = 1 if ARGS.FORWARD_SELECT else ARGS.NUM_FEATURES
    forward_select = None
    results = None
    for num_features in range(forward_initval, ARGS.NUM_FEATURES + 1):

        if sim is None:
            ylabeler = Labeler(
                seqrecord_get_ic50s,
                lambda row: is_refseq(row) or False, # TODO: again filtration function
                lambda x: x > ARGS.IC50,
                ARGS.AUTOBALANCE
            )
            y, ic50 = ylabeler(alignment)
            assert y.shape[0] == x.shape[0], "number of classes doesn't match the data: %d vs %d" % (y.shape[0], x.shape[0])
            assert(
                (ic50 is None and not ARGS.AUTOBALANCE) or
                (ic50 is not None and ARGS.AUTOBALANCE)
            )
            if ARGS.AUTOBALANCE:
                ARGS.IC50 = ic50

        # simulations, ho!
        for i in range(sim.runs if sim is not None else 1):

            # here is where the sequences must be generated for the random sequence and random epitope simulations
            if sim is not None:
                alignment = sim.generate_sequences(
                    N=ARGS.SIM_EPI_N,
                    idfmt='%s|||',
                    noise=ARGS.SIM_EPI_NOISE,
                    mutation_rate=ARGS.SIM_EPI_MUT_RATE,
                    alphabet=alph
                )

                builder = DataBuilder(
                    alignment,
                    alph,
                    refidx,
                    filter
                )
                x = builder(alignment, refidx)
                colnames = builder.labels

                # simulates the epitope and assigns the appropriate class
                epi_def = sim.simulate_epitope(
                    alignment,
                    alph,
                    colnames,
                    ARGS.SIM_EPI_SIZE,
                    ARGS.SIM_EPI_PERCENTILE,
                )

                if epi_def is not None:
                    print('********************* SIMULATED EPITOPE DESCRIPTION (%d) *********************\n' % ARGS.SIM_EPI_SIZE, file=sys.stdout)
                    print('%s\n' % str(epi_def), file=sys.stdout)

            if ARGS.MRMR_NORMALIZE:
                DiscreteMrmr._NORMALIZED = True

            crossvalidator = SelectingNestedCrossValidator(
                classifier_cls=LinearSvm,
                selector_cls=DiscreteMrmr,
                folds=ARGS.CV_FOLDS,
                gridsearch_kwargs={ 'C': C_range(*ARGS.LOG2C) },
                classifier_kwargs={ 'bias': True },
                selector_kwargs={
                    'num_features': num_features,
                    'method': ARGS.MRMR_METHOD,
                    'threshold': ARGS.SIMILAR
                },
                validator_cls=CrossValidator,
                validator_kwargs={
                    'folds': ARGS.CV_FOLDS-1,
                    'scorer_cls': DiscretePerfStats,
                    'scorer_kwargs': { 'optstat': ARGS.OPTSTAT }
                },
                scorer_cls=DiscretePerfStats,
                scorer_kwargs={ 'optstat': ARGS.OPTSTAT },
                weights_func='weights' # we MUST specify this or it will be set to lambda: None
            )

            new_results = crossvalidator.crossvalidate(
                x, y,
                classifier_kwargs={},
                extra=extract_feature_weights_similar if similar else extract_feature_weights
            )

            if results is not None and new_results.stats.get() <= results.stats.get():
                break

            results = new_results
            forward_select = num_features

    # the alignment reflects the number of sequences either naturally,
    # or through SIM_EPI_N, which reflects the natural number anyway, less the refseq
    ret = IdepiResults(similar)
    ret.cv_results(results, colnames)
    ret.metadata(ARGS, len(alignment)-1, np.mean(y), antibody, forward_select)

    print(ret.dumps(), file=ARGS.OUTPUT)

    finalize_args(ARGS)

    return ret


if __name__ == '__main__':
    sys.exit(main())
