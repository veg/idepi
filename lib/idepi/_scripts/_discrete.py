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
from operator import itemgetter
from os import close, remove
from os.path import abspath, join, split
from random import seed
from tempfile import mkstemp
from warnings import warn

import numpy as np

from Bio import AlignIO

from BioExt import hxb2

from idepi import (
    Alphabet,
    ClassExtractor,
    DumbSimulation,
    LinearSvm,
    MarkovSimulation,
    NaiveFilter,
    Simulation,
    alignment_identify_ref,
    cv_results_to_output,
    extract_feature_weights,
    extract_feature_weights_similar,
    fasta_json_desc,
    fix_hxb2_seq,
    generate_alignment,
    input_data,
    is_refseq,
    make_output_meta,
    pretty_fmt_results,
    seqrecord_get_ic50s,
    set_util_params,
    __file__ as _idepi_file,
    __version__ as _idepi_version
)

from mrmr import DiscreteMrmr

from pyxval import CrossValidator, DiscretePerfStats, SelectingNestedCrossValidator


__VERSION__ = _idepi_version

_IDEPI_PATH       = split(abspath(_idepi_file))[0]
_HXB2_DNA_FASTA   = join(_IDEPI_PATH, 'data', 'hxb2_dna.fa')
_HXB2_AMINO_FASTA = join(_IDEPI_PATH, 'data', 'hxb2_pep.fa')

_DEFAULT_NUM_FEATURES = 10

_RAND_SEQ_STOCKHOLM = join(_IDEPI_PATH, 'data', 'randhivenvpep_final.sto')

# strip the _TEST variables because of the beginning and trailing newlines
_TEST_DNA_STO = '''# STOCKHOLM 1.0
1||A|1       AUGAUUCCCGACUUUAAANNN
2||A|21      AUGAUUCCCGACUUUAAANNNCAC
3||A|50      AUGAUUCCCAAANNNCAC
4||B|0.5     AUGCCCGACUUUAAACAC
HXB2_env     AUGCCCGACUUUAAACAC
//'''.strip()

_TEST_AMINO_STO = '''# STOCKHOLM 1.0
1||A|1        MIPDFKX-
2||A|21       MIPDFKXH
3||A|50       MIP--KXH
4||B|0.5      .MPDFKH-
HXB2_env      -MPDFKH-
//'''.strip()

_TEST_AMINO_NAMES = ['0aM', '0a[]', 'M1I', 'M1M', 'P2P', 'D3D', 'D3[]', 'F4F', 'F4[]', 'K5K', 'H6H', 'H6X', '6aH', '6a[]']
_TEST_STANFEL_NAMES = ['0a[ACGILMPSTV]', '0a[]', 'M1[ACGILMPSTV]', 'P2[ACGILMPSTV]', 'D3[DENQ]', 'D3[]', \
                       'F4[FWY]', 'F4[]', 'K5[HKR]', 'H6[HKR]', 'H6[X]', '6a[HKR]', '6a[]']

_TEST_Y = np.array([1,0,0,1])

_TEST_AMINO_X = np.array([[1,0,1,0,1,1,0,1,0,1,0,1,0,1],
                          [1,0,1,0,1,1,0,1,0,1,0,1,1,0],
                          [1,0,1,0,1,0,1,0,1,1,0,1,1,0],
                          [0,1,0,1,1,1,0,1,0,1,1,0,0,1]])

_TEST_STANFEL_X = np.array([[1,0,1,1,1,0,1,0,1,0,1,0,1],
                            [1,0,1,1,1,0,1,0,1,0,1,1,0],
                            [1,0,1,1,0,1,0,1,1,0,1,1,0],
                            [0,1,1,1,1,0,1,0,1,1,0,0,1]])

ARGS = None


def setup_option_parser():

    parser = OptionParser(usage = '%prog [options] ANTIBODY')

    #                 option             action = 'store'      callback                  type           nargs=1  dest
    parser.add_option('--hmmalign',                                                      type='string',          dest='HMMER_ALIGN_BIN')
    parser.add_option('--hmmbuild',                                                      type='string',          dest='HMMER_BUILD_BIN')
    parser.add_option('--hmmiter',                                                       type='int',             dest='HMMER_ITER')
    parser.add_option('--log',           action='callback',    callback=optparse_log,    type='string',          dest='LOGGING')
    parser.add_option('--mrmr',          action='store_false',                                                   dest='MAXREL')
    parser.add_option('--mrmrmethod',                                                    type='string',          dest='MRMR_METHOD')
    parser.add_option('--maxrel',        action='store_true',                                                    dest='MAXREL')
    parser.add_option('--normalizemrmr', action='store_true',                                                    dest='MRMR_NORMALIZE')
    parser.add_option('--filter',        action='callback',    callback=optparse_csv,    type='string',          dest='FILTER')
    parser.add_option('--clonal',        action='store_true',                                                    dest='CLONAL')
    parser.add_option('--numfeats',                                                      type='int',             dest='NUM_FEATURES')
    parser.add_option('--forward',       action='store_true',                                                    dest='FORWARD_SELECT')
    parser.add_option('--subtypes',      action='callback',    callback=optparse_csv,    type='string',          dest='SUBTYPES')
    parser.add_option('--log2c',         action='callback',    callback=optparse_csv,    type='string',          dest='LOG2C')
    parser.add_option('--weighting',     action='store_true',                                                    dest='WEIGHTING')
    parser.add_option('--accuracy',      action='store_true',                                                    dest='ACCURACY')
    parser.add_option('--ppv',           action='store_true',                                                    dest='PPV')
    parser.add_option('--precision',     action='store_true',                                                    dest='PPV')
    parser.add_option('--npv'     ,      action='store_true',                                                    dest='NPV')
    parser.add_option('--sensitivity',   action='store_true',                                                    dest='SENSITIVITY')
    parser.add_option('--recall',        action='store_true',                                                    dest='SENSITIVITY')
    parser.add_option('--specificity',   action='store_true',                                                    dest='SPECIFICITY')
    parser.add_option('--tnr',           action='store_true',                                                    dest='SPECIFICITY')
    parser.add_option('--fscore',        action='store_true',                                                    dest='FSCORE')
    parser.add_option('--amino',         action='store_true',                                                    dest='AMINO')
    parser.add_option('--dna',           action='store_true',                                                    dest='DNA')
    parser.add_option('--stanfel',       action='store_true',                                                    dest='STANFEL')
    parser.add_option('--cv',                                                            type='int',             dest='CV_FOLDS')
    parser.add_option('--loocv',         action='store_true',                                                    dest='LOOCV')
    parser.add_option('--targets',       action='callback',    callback=optparse_csv,    type='string',          dest='TARGETS')
    parser.add_option('--maxcon',                                                        type='float',           dest='MAX_CONSERVATION')
    parser.add_option('--maxgap',                                                        type='float',           dest='MAX_GAP_RATIO')
    parser.add_option('--mincon',                                                        type='float',           dest='MIN_CONSERVATION')
    parser.add_option('--ic50gt',                                                        type='float',           dest='IC50GT')
    parser.add_option('--ic50lt',                                                        type='float',           dest='IC50LT')
    parser.add_option('--data',          action='callback',    callback=optparse_data,   type='string',          dest='DATA')
    parser.add_option('--refseq',                                                        type='string',          dest='REFSEQ')
    parser.add_option('--ids',           action='callback',    callback=optparse_csv,    type='string',          dest='HXB2_IDS')
    parser.add_option('--test',          action='store_true',                                                    dest='TEST')
    parser.add_option('--sim',                                                           type='string',          dest='SIM')
    parser.add_option('--simruns',                                                       type='int',             dest='SIM_RUNS')
    parser.add_option('--simepisize',                                                    type='int',             dest='SIM_EPI_SIZE')
    parser.add_option('--simepimutrate',                                                 type='float',           dest='SIM_EPI_MUT_RATE')
    parser.add_option('--simepiseqnum',                                                  type='int',             dest='SIM_EPI_N')
    parser.add_option('--simepinoise',                                                   type='float',           dest='SIM_EPI_NOISE')
    parser.add_option('--simepiperc',                                                    type='float',           dest='SIM_EPI_PERCENTILE')
    parser.add_option('--seed',                                                          type='int',             dest='RAND_SEED')
    parser.add_option('--phylofilt',     action='store_true',                                                    dest='PHYLOFILTER')
    parser.add_option('-o', '--output',                                                  type='string',          dest='OUTPUT')

    parser.set_defaults(LOGGING            = None)
    parser.set_defaults(HMMER_ALIGN_BIN    = 'hmmalign')
    parser.set_defaults(HMMER_BUILD_BIN    = 'hmmbuild')
    parser.set_defaults(HMMER_ITER         = 8)
    parser.set_defaults(MRMR_METHOD        = 'MID')
    parser.set_defaults(MRMR_NORMALIZE     = False)
    parser.set_defaults(MAXREL             = False)
    parser.set_defaults(FILTER             = [])
    parser.set_defaults(CLONAL             = False)
    parser.set_defaults(NUM_FEATURES       = -1)
    parser.set_defaults(FORWARD_SELECT     = False)
    parser.set_defaults(SUBTYPES           = [])
    parser.set_defaults(WEIGHTING          = False)
    parser.set_defaults(LOG2C              = ['-5', '15', '0.25'])
    parser.set_defaults(ACCURACY           = False)
    parser.set_defaults(PPV                = False)
    parser.set_defaults(NPV                = False)
    parser.set_defaults(SENSITIVITY        = False)
    parser.set_defaults(SPECIFICITY        = False)
    parser.set_defaults(FSCORE             = False)
    parser.set_defaults(AMINO              = False)
    parser.set_defaults(DNA                = False)
    parser.set_defaults(STANFEL            = False)
    parser.set_defaults(CV_FOLDS           = 5)
    parser.set_defaults(LOOCV              = False)
    parser.set_defaults(TARGETS            = ['gt'])
    parser.set_defaults(MAX_CONSERVATION   = 1. ) # 93.)
    parser.set_defaults(MAX_GAP_RATIO      = 0.1 ) # 93.)
    parser.set_defaults(MIN_CONSERVATION   = 1. ) # 33.)
    parser.set_defaults(IC50GT             = 20.)
    parser.set_defaults(IC50LT             = 2.)
    parser.set_defaults(DATA               = input_data(join(_IDEPI_PATH, 'data', 'allneuts.sqlite3')))
    parser.set_defaults(REFSEQ             = hxb2.env.load())
    parser.set_defaults(HXB2_IDS           = ['HXB2_env'])
    parser.set_defaults(SIM                = '') # can be 'randtarget' for now
    parser.set_defaults(SIM_RUNS           = 1) # can be 'randtarget' for now
    parser.set_defaults(SIM_EPI_SIZE       = 10)
    parser.set_defaults(SIM_EPI_MUT_RATE   = 0.01)
    parser.set_defaults(SIM_EPI_N          = None) # default is to use len(seqrecords)
    parser.set_defaults(SIM_EPI_NOISE      = 0.08)
    parser.set_defaults(SIM_EPI_PERCENTILE = 0.5)
    parser.set_defaults(RAND_SEED          = 42) # magic number for determinism
    parser.set_defaults(PHYLOFILTER        = False)
    parser.set_defaults(OUTPUT             = None)

    return parser


def run_tests():
    # set these to this so we don't exclude anything (just testing file generation and parsing)
    ARGS.NUM_FEATURES = 15 # should be enough, the number is known to be 13
    ARGS.MRMR_METHOD = DiscreteMrmr.MID
    ARGS.MAX_CONSERVATION = 1.0
    ARGS.MAX_GAP_RATIO    = 1.0
    ARGS.MIN_CONSERVATION = 1.0
    ARGS.IC50 = 20.

    # if we don't do this, DOOMBUNNIES
    set_util_params(ARGS.REFSEQ_IDS, ARGS.IC50)

    fd, sto_filename = mkstemp(); close(fd)

    try:
        fh = open(sto_filename, 'w')
        print(_TEST_AMINO_STO, file=fh)
        fh.close()

        alignment = AlignIO.read(sto_filename, 'stockholm')
        refidx = alignment_identify_ref(alignment, is_refseq)

        for ARGS.ALPHABET in (Alphabet.AMINO, Alphabet.STANFEL):

            if ARGS.ALPHABET == Alphabet.STANFEL:
                _TEST_NAMES = _TEST_STANFEL_NAMES
                _TEST_X = _TEST_STANFEL_X
            else:
                _TEST_NAMES = _TEST_AMINO_NAMES
                _TEST_X = _TEST_AMINO_X

            alph = Alphabet(mode=ARGS.ALPHABET)

            colfilter = NaiveFilter(
                alph,
                ARGS.MAX_CONSERVATION,
                ARGS.MIN_CONSERVATION,
                ARGS.MAX_GAP_RATIO,
                refidx=refidx,
                skip_func=lambda x: False # TODO: add the appropriate filter function based on the args here
            )
            colnames, x = colfilter.learn(alignment, {})

            # test the feature names portion
            try:
                assert(len(colnames) == len(_TEST_NAMES))
            except AssertionError:
                raise AssertionError('gen:   %s\ntruth: %s' % (colnames, _TEST_NAMES))

            for name in _TEST_NAMES:
                try:
                    assert(name in colnames)
                except AssertionError:
                    raise AssertionError('ERROR: \'%s\' not found in %s' % (name, ', '.join(colnames)))

            assert(np.all(_TEST_X == x))

            # test mRMR and LSVM file generation
            yextractor = ClassExtractor(
                seqrecord_get_ic50s,
                lambda row: is_refseq(row) or False, # TODO: again filtration function
                lambda x: x >  ARGS.IC50,
                False
            )
            y, ic50 = yextractor.extract(alignment)

            assert(np.all(_TEST_Y == y))

            if ARGS.MRMR_NORMALIZE:
                DiscreteMrmr._NORMALIZED = True

            # generate and test the mRMR portion
            mrmr = DiscreteMrmr(
                num_features=ARGS.NUM_FEATURES,
                method=ARGS.MRMR_METHOD
            )

            mrmr.select(x, y)

            x = mrmr.subset(x)

            lsvm = LinearSvm()
            lsvm.learn(x, y)

    finally:
        remove(sto_filename)

    print('ALL TESTS PASS', file=sys.stderr)


def main(argv=sys.argv):
    np.seterr(all='raise')

    global ARGS

    # so some option parsing
    parser = init_args("XXX: description")
    ARGS = parser.parse_args(argv)

    finalize_args(ARGS)

    # do some argument parsing
    if ARGS.TEST:
        run_tests()
        return 0

    if ARGS.RAND_SEED is not None:
        seed(ARGS.RAND_SEED)

    similar = True
    if ARGS.MAXREL:
        similar = False

    # XXX autobalance option directly
    autobalance = False
    if ARGS.TARGETS == ['auto']:
        ARGS.TARGETS = ['ge']
        autobalance = True

    # validate the subtype option
    valid_subtypes = sorted(ARGS.DATA.subtypes, key=lambda x: x.strip().upper())
    for subtype in ARGS.SUBTYPES:
        if subtype not in valid_subtypes:
            msg = "'%s' not in the list of permitted subtypes: %s" % (
                subtype,
                ', '.join("'%s'" % st.strip() for st in valid_subtypes)
            )
            raise ArgumentTypeError(msg)

    # validate the antibody argument, currently a hack exists to make PG9/PG16 work
    # TODO: Fix pg9/16 hax
    antibody = args.antibody.strip()
    valid_antibodies = sorted(ARGS.DATA.antibodies, key=lambda x: x.strip())
    if antibody not in valid_antibodies:
        if ' ' + antibody not in valid_antibodies:
            msg = "'%s' not in the list of permitted antibodies: %s" % (
                antibody,
                ', '.join("'%s'" % ab.strip() for ab in valid_antibodies)
            )
            raise ArgumentTypeError(msg)
        else:
            antibody = ' ' + antibody

    if ARGS.SIM in (Simulation.EPITOPE, Simulation.SEQUENCE) and ARGS.DNA:
        raise ArgumentTypeError('randseq simulation target not compatible with DNA alphabet')

    if len(ARGS.FILTER) != 0:
        if ARGS.NUM_FEATURES > len(ARGS.FILTER):
            ARGS.NUM_FEATURES = len(ARGS.FILTER)
            warn('clamping --numfeats to sizeof(--filter) = %d' % ARGS.NUM_FEATURES)

    # convert the hxb2 reference to amino acid, including loop definitions
    fix_hxb2_seq(ARGS)

    # set the util params
    set_util_params(ARGS.REFSEQ_IDS, ARGS.IC50)

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=ARGS.ALPHABET)

    sim = None
    if ARGS.SIM is not None:
        if ARGS.SIM == Simulation.DUMB:
            sim = DumbSimulation(ARGS.SIM_RUNS, Simulation.EPITOPE, str(ARGS.REFSEQ.seq))
        elif ARGS.SIM in (Simulation.EPITOPE, Simulation.SEQUENCE):
            sim = MarkovSimulation(ARGS.SIM_RUNS, ARGS.SIM, _RAND_SEQ_STOCKHOLM)
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
        __VERSION__
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

    alignment, refseq_offs = generate_alignment(seqrecords, alignment_basename, is_refseq, ARGS)
    refidx = alignment_identify_ref(alignment, is_refseq)
    colfilter = NaiveFilter(
        alph,
        ARGS.MAX_CONSERVATION,
        ARGS.MIN_CONSERVATION,
        ARGS.MAX_GAP_RATIO,
        refidx=refidx,
        skip_func=lambda x: False, # TODO: add the appropriate filter function based on the args here
        loop_defs=sorted(fasta_json_desc(ARGS.REFSEQ).get('loops', {}).items(), key=itemgetter(0))
    )

    if sim is None:
        colnames, x = colfilter.learn(alignment, refseq_offs)
        # TODO: I don't think we need to binarize the colnames here, though we can if we want.
        # I need to think more about how to properly handle this case.
    else:
        if ARGS.SIM_EPI_N is None:
            ARGS.SIM_EPI_N = len(seqrecords)

    # compute features
    forward_initval = 1 if ARGS.FORWARD_SELECT else ARGS.NUM_FEATURES
    forward_select = None
    results = None
    for num_features in range(forward_initval, ARGS.NUM_FEATURES + 1):

        if sim is None:
            yextractor = ClassExtractor(
                seqrecord_get_ic50s,
                lambda row: is_refseq(row) or False, # TODO: again filtration function
                lambda x: x > ARGS.IC50,
                ARGS.AUTOBALANCE
            )
            y, ic50 = yextractor.extract(alignment)
            assert y.shape[0] == x.shape[0], "number of classes doesn't match the data: %d vs %d" % (y.shape[0], x.shape[0])
            assert(
                (ic50 is None and not ARGS.AUTOBALANCE) or
                (ic50 is not None and ARGS.AUTOBALANCE)
            )
            if autobalance:
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

                colnames, x = colfilter.learn(alignment, {}) # XXX: refseq_offs needed here?

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

            C_begin, C_end, C_step = ARGS.LOG2C
            recip = 1
            if isinstance(C_step, float):
                recip = 1. / C_step
                C_begin, C_end = int(recip * C_begin), int(recip * C_end)
                C_step = 1
            C_range = [pow(2., float(C) / recip) for C in range(C_begin, C_end + 1, C_step)]

            if ARGS.MRMR_NORMALIZE:
                DiscreteMrmr._NORMALIZED = True

            crossvalidator = SelectingNestedCrossValidator(
                classifier_cls=LinearSvm,
                selector_cls=DiscreteMrmr,
                folds=ARGS.CV_FOLDS,
                gridsearch_kwargs={ 'C': C_range },
                classifier_kwargs={ 'bias': True },
                selector_kwargs={
                    'num_features': num_features,
                    'method': ARGS.MRMR_METHOD
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
    meta = make_output_meta(ARGS, len(alignment)-1, np.mean(y), antibody, forward_select)
    ret = cv_results_to_output(results, colnames, meta, similar)

    if isinstance(ARGS.OUTPUT, str):
        with open(ARGS.OUTPUT, 'w') as fh:
            print(pretty_fmt_results(ret, similar), file=fh)
    else:
        print(pretty_fmt_results(ret, similar))

    return 0


if __name__ == '__main__':
    sys.exit(main())
