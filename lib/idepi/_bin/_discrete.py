#!/usr/bin/env python2.7
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

import logging, sqlite3, sys

from math import ceil, copysign, log10, sqrt
from operator import itemgetter
from optparse import OptionParser
from os import close, remove, rename
from os.path import abspath, basename, exists, join, split, splitext
from random import gauss, random, seed
from re import sub, match
from tempfile import mkstemp

import numpy as np

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from idepi import (Alphabet, ClassExtractor, DumbSimulation, LinearSvm, MarkovSimulation, NaiveFilter,
                   PhyloFilter, SeqTable, Simulation, collect_seqrecords_from_db, cv_results_to_output,
                   extract_feature_weights, generate_alignment, get_valid_antibodies_from_db, IDEPI_LOGGER,
                   get_valid_subtypes_from_db, is_HXB2, make_output_meta, pretty_fmt_results,
                   seqrecord_to_ic50s, set_util_params, __file__ as _idepi_file,
                   __version__ as _idepi_version)

from mrmr import MRMR_LOGGER, DiscreteMrmr, PhyloMrmr

from pyxval import PYXVAL_LOGGER, CrossValidator, DiscretePerfStats, SelectingNestedCrossValidator


__VERSION__ = _idepi_version

_IDEPI_PATH       = split(abspath(_idepi_file))[0]
_HXB2_DNA_FASTA   = join(_IDEPI_PATH, 'data', 'hxb2_dna.fa')
_HXB2_AMINO_FASTA = join(_IDEPI_PATH, 'data', 'hxb2_pep.fa')

_DEFAULT_NUM_FEATURES = 10

_RAND_SEQ_STOCKHOLM = join(_IDEPI_PATH, 'data', 'randhivenvpep_final.sto')

_PHYLOFILTER_BATCHFILE = join(_IDEPI_PATH, 'hyphy', 'CorrectForPhylogeny.bf')

# strip the _TEST variables because of the beginning and trailing newlines
_TEST_DNA_STO = '''# STOCKHOLM 1.0
1||A|1       AUGAUUCCCGACUUUAAANNN
2||A|21      AUGAUUCCCGACUUUAAANNNCAC
3||A|50      AUGAUUCCCAAANNNCAC
4||B|0.5     AUGCCCGACUUUAAACAC
gi|9629357   AUGCCCGACUUUAAACAC
//'''.strip()

_TEST_AMINO_STO = '''# STOCKHOLM 1.0
1||A|1        MIPDFKX-
2||A|21       MIPDFKXH
3||A|50       MIP--KXH
4||B|0.5      .MPDFKH-
gi|9629363    -MPDFKH-
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

OPTIONS = None


def optparse_extend(option, opt_str, value, parser):
    if getattr(parser.values, option.dest, None) is None:
        setattr(parser.values, option.dest, [])
    getattr(parser.values, option.dest).extend(value)


def optparse_csv(option, opt_str, value, parser):
    setattr(parser.values, option.dest, value.split(','))


def setup_option_parser():

    parser = OptionParser(usage = '%prog [options] ANTIBODY')

    #                 option             action = 'store'      callback                  type           nargs=1  dest
    parser.add_option('--hmmalign',                                                      type='string',          dest='HMMER_ALIGN_BIN')
    parser.add_option('--hmmbuild',                                                      type='string',          dest='HMMER_BUILD_BIN')
    parser.add_option('--hmmiter',                                                       type='int',             dest='HMMER_ITER')
    parser.add_option('--log',           action='store_true',                                                    dest='LOGGING')
    parser.add_option('--mrmr',          action='store_false',                                                   dest='MAXREL')
    parser.add_option('--mrmrmethod',                                                    type='string',          dest='MRMR_METHOD')
    parser.add_option('--maxrel',        action='store_true',                                                    dest='MAXREL')
    parser.add_option('--normalizemrmr', action='store_true',                                                    dest='MRMR_NORMALIZE')
    parser.add_option('--filter',        action='callback',    callback=optparse_csv,    type='string',          dest='FILTER')
    parser.add_option('--numfeats',                                                      type='int',             dest='NUM_FEATURES')
    parser.add_option('--forward',       action='store_true',                                                    dest='FORWARD_SELECT')
    parser.add_option('--subtypes',      action='callback',    callback=optparse_csv,    type='string',          dest='SUBTYPES')
    parser.add_option('--nested',        action='store_true',                                                    dest='NESTED')
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
    parser.add_option('--neuts',                                                         type='string',          dest='NEUT_SQLITE3_DB')
    parser.add_option('--refseq',                                                        type='string',          dest='REFSEQ_FASTA')
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
    parser.add_option('--output',                                                        type='string',          dest='OUTPUT')

    parser.set_defaults(LOGGING            = False)
    parser.set_defaults(HMMER_ALIGN_BIN    = 'hmmalign')
    parser.set_defaults(HMMER_BUILD_BIN    = 'hmmbuild')
    parser.set_defaults(HMMER_ITER         = 8)
    parser.set_defaults(MRMR_METHOD        = 'MID')
    parser.set_defaults(MRMR_NORMALIZE     = False)
    parser.set_defaults(MAXREL             = False)
    parser.set_defaults(FILTER             = [])
    parser.set_defaults(NUM_FEATURES       = -1)
    parser.set_defaults(FORWARD_SELECT     = False)
    parser.set_defaults(SUBTYPES           = [])
    parser.set_defaults(NESTED             = True)
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
    parser.set_defaults(NEUT_SQLITE3_DB    = join(_IDEPI_PATH, 'data', 'allneuts.sqlite3'))
    parser.set_defaults(REFSEQ_FASTA       = _HXB2_AMINO_FASTA)
    parser.set_defaults(HXB2_IDS           = ['9629357', '9629363'])
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
    OPTIONS.NUM_FEATURES = 15 # should be enough, the number is known to be 13
    OPTIONS.MAXREL = False
    OPTIONS.DNA = False
    OPTIONS.MAX_CONSERVATION = 1.0
    OPTIONS.MAX_GAP_RATIO    = 1.0
    OPTIONS.MIN_CONSERVATION = 1.0
    OPTIONS.IC50GT = 20.
    OPTIONS.IC50LT = 2.
    OPTIONS.TARGETS = ['lt']

    # if we don't do this, DOOMBUNNIES
    set_util_params(OPTIONS.HXB2_IDS, OPTIONS.IC50GT, OPTIONS.IC50LT)

    fd, sto_filename = mkstemp(); close(fd)

    try:
        fh = open(sto_filename, 'w')
        print(_TEST_AMINO_STO, file=fh)
        fh.close()

        alignment = AlignIO.read(sto_filename, 'stockholm')

        for OPTIONS.STANFEL in (True, False):

            if OPTIONS.STANFEL:
                OPTIONS.AMINO = False
                _TEST_NAMES = _TEST_STANFEL_NAMES
                _TEST_X = _TEST_STANFEL_X
            else:
                OPTIONS.AMINO = True
                _TEST_NAMES = _TEST_AMINO_NAMES
                _TEST_X = _TEST_AMINO_X

            alph = Alphabet(Alphabet.STANFEL if OPTIONS.STANFEL else Alphabet.DNA if OPTIONS.DNA else Alphabet.AMINO)

            colfilter = NaiveFilter(
                alph,
                OPTIONS.MAX_CONSERVATION,
                OPTIONS.MIN_CONSERVATION,
                OPTIONS.MAX_GAP_RATIO,
                is_HXB2,
                lambda x: False # TODO: add the appropriate filter function based on the args here
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
            for target in OPTIONS.TARGETS:
                yextractor = ClassExtractor(
                    seqrecord_to_ic50s,
                    lambda row: is_HXB2(row) or False, # TODO: again filtration function
                    lambda x: x < OPTIONS.IC50LT if target == 'lt' else x > OPTIONS.IC50GT
                )
                y = yextractor.extract(alignment)

                assert(np.all(_TEST_Y == y))

                if OPTIONS.MRMR_NORMALIZE:
                    DiscreteMrmr._NORMALIZED = True

                # generate and test the mRMR portion
                mrmr = DiscreteMrmr(
                    num_features=OPTIONS.NUM_FEATURES,
                    method=DiscreteMrmr.MAXREL if OPTIONS.MAXREL \
                      else DiscreteMrmr.MID if OPTIONS.MRMR_METHOD == 'MID' \
                      else DiscreteMrmr.MIQ
                )

                mrmr.select(x, y)

                x = mrmr.subset(x)

                lsvm = LinearSvm()
                lsvm.learn(x, y)

    finally:
        remove(sto_filename)

    print('ALL TESTS PASS', file=sys.stderr)


def fix_hxb2_fasta():
    '''If DNA mode was selected but the AMINO reference sequence is still in place, fix it'''
    if OPTIONS.DNA == True and OPTIONS.REFSEQ_FASTA == _HXB2_AMINO_FASTA:
        OPTIONS.REFSEQ_FASTA = _HXB2_DNA_FASTA


def main(argv=sys.argv):
    np.seterr(all='raise')

    global OPTIONS

    # so some option parsing
    option_parser = setup_option_parser()
    (OPTIONS, args) = option_parser.parse_args(argv)

    if OPTIONS.LOGGING:
        logging.getLogger(IDEPI_LOGGER).setLevel(logging.DEBUG)
        logging.getLogger(MRMR_LOGGER).setLevel(logging.DEBUG)
        logging.getLogger(PYXVAL_LOGGER).setLevel(logging.DEBUG)
        try:
            from fakemp import FAKEMP_LOGGER
            logging.getLogger(FAKEMP_LOGGER).setLevel(logging.DEBUG)
        except ImportError:
            pass

    # do some argument parsing
    if OPTIONS.TEST:
        run_tests()
        return 0

    if OPTIONS.SIM != '':
        if OPTIONS.SIM == 'randdumbepi':
            OPTIONS.SIM = Simulation.DUMB
        elif OPTIONS.SIM == 'randepi':
            OPTIONS.SIM == Simulation.EPITOPE
        elif OPTIONS.SIM == 'randseq':
            OPTIONS.SIM == Simulation.SEQUENCE
        elif OPTIONS.SIM == 'randtarget':
            OPTIONS.SIM == Simulation.TARGET
        else:
            option_parser.error('option --sim takes one of %s' % ', '.join(Simulation.VALUES))

    if OPTIONS.RAND_SEED is not None:
        seed(OPTIONS.RAND_SEED)

    if len(args) != 2:
        option_parser.error('ANTIBODY is a required argument')

    if not set(OPTIONS.TARGETS).issubset(set(['lt', 'gt'])):
        option_parser.error('option --targets takes either or both: lt gt')

    if not OPTIONS.MRMR_METHOD in ('MIQ', 'MID'):
        option_parser.error('option --mrmrmethod takes either MIQ or MID')

    # check to make sure our mode is exclusive, and set the default (AMINO) if none is set
    alphoptlen = sum(1 for v in (OPTIONS.AMINO, OPTIONS.DNA, OPTIONS.STANFEL) if v)
    if alphoptlen > 1:
        option_parser.error('options --amino, --dna, and --stanfel are mutually exclusive')
    elif alphoptlen == 0:
        OPTIONS.AMINO = True

    if sum(1 for v in (OPTIONS.ACCURACY, OPTIONS.PPV, OPTIONS.NPV, OPTIONS.SENSITIVITY, OPTIONS.SPECIFICITY, OPTIONS.FSCORE) if v) > 1:
        option_parser.error('options --accuracy, --ppv/--precision, --npv, --sensitivity/--recall, --specificity/--tnr, --fscore are mutually exclusive')

    try:
        if len(OPTIONS.LOG2C) != 3 or float(OPTIONS.LOG2C[2]) <= 0.:
            raise ValueError
        OPTIONS.LOG2C = [int(OPTIONS.LOG2C[0]), int(OPTIONS.LOG2C[1]), float(OPTIONS.LOG2C[2])]
    except ValueError as e:
        option_parser.error('option --log2c takes an argument of the form C_BEGIN,C_END,C_STEP where C_STEP must be > 0')

    # validate the antibody argument, currently a hack exists to make PG9/PG16 work
    # TODO: Fix pg9/16 hax
    antibody = args[1]
    valid_antibodies = sorted(get_valid_antibodies_from_db(OPTIONS.NEUT_SQLITE3_DB), key = lambda x: x.strip())
    if antibody not in valid_antibodies:
        if ' ' + antibody not in valid_antibodies:
            option_parser.error('%s not in the list of permitted antibodies: \n  %s' % (antibody, '\n  '.join(ab.strip() for ab in valid_antibodies)))
        else:
            antibody = ' ' + antibody

    # validate the subtype option
    valid_subtypes = sorted(get_valid_subtypes_from_db(OPTIONS.NEUT_SQLITE3_DB), key = lambda x: x.strip().upper())
    for subtype in OPTIONS.SUBTYPES:
        if subtype not in valid_subtypes:
            option_parser.error('%s not in the list of permitted subtypes: \n  %s' % (subtype, '\n  '.join(st.strip() for st in valid_subtypes)))

    if OPTIONS.SIM in (Simulation.EPITOPE, Simulation.SEQUENCE) and OPTIONS.DNA:
        option_parser.error('randseq simulation target not compatible with DNA mode')

    if OPTIONS.SIM_EPI_MUT_RATE < 0. or OPTIONS.SIM_EPI_MUT_RATE > 1.:
        option_parser.error('--simepimutrate must be between 0.0 and 1.0')

    if abs(OPTIONS.SIM_EPI_NOISE) > 1.:
        option_parser.error('--simepinoise shouldn\'t really ever be greater than 1.0')

    if OPTIONS.SIM_EPI_N is not None and OPTIONS.SIM_EPI_N < 1:
        option_parser.error('--simepiseqnum must be greater than 0')

    if OPTIONS.SIM_EPI_PERCENTILE < 0. or OPTIONS.SIM_EPI_PERCENTILE > 1.:
        option_parser.error('--simepiperc must be betweeen 0.0 and 1.0')

    if len(OPTIONS.FILTER) != 0:
        if OPTIONS.NUM_FEATURES > len(OPTIONS.FILTER):
            OPTIONS.NUM_FEATURES = len(OPTIONS.FILTER)
            print('warning: clamping --numfeats to sizeof(--filter) = %d' % OPTIONS.NUM_FEATURES, file=sys.stderr)
    else: # len(OPTIONS.FILTER) == 0
        if OPTIONS.NUM_FEATURES == -1:
            OPTIONS.NUM_FEATURES = _DEFAULT_NUM_FEATURES

    if OPTIONS.MRMR_NORMALIZE and OPTIONS.PHYLOFILTER:
        print("mRMR normalization and phylofiltering are incompatible, disabling mRMR normalization", file=sys.stderr)
        OPTIONS.MRMR_NORMALIZE = False

    # destroy the parser because optparse docs recommend it
    option_parser.destroy()

    # use the default DNA HXB2 Reference seq if we define --dna but don't give a new default HXB2 Reference seq
    fix_hxb2_fasta()

    # set the util params
    set_util_params(OPTIONS.HXB2_IDS, OPTIONS.IC50GT, OPTIONS.IC50LT)

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=Alphabet.STANFEL if OPTIONS.STANFEL else Alphabet.DNA if OPTIONS.DNA else Alphabet.AMINO)

    sim = None
    if OPTIONS.SIM != '':
        if OPTIONS.SIM == Simulation.DUMB:
            hxb2fh = open(OPTIONS.REFSEQ_FASTA)
            for record in SeqIO.parse(hxb2fh, 'fasta'):
                hxb2seq = str(record.seq)
                break
            hxb2fh.close()
            sim = DumbSimulation(OPTIONS.SIM_RUNS, Simulation.EPITOPE, hxb2seq)
        elif OPTIONS.SIM in (Simulation.EPITOPE, Simulation.SEQUENCE):
            sim = MarkovSimulation(OPTIONS.SIM_RUNS, OPTIONS.SIM, _RAND_SEQ_STOCKHOLM)
        else:
            raise ValueError('Unknown simulation type `%s\'' % OPTIONS.SIM)

    ab_basename = '%s%s_%s' % (antibody, '_randseq' if sim is not None and sim.mode == Simulation.SEQUENCE else '', 'dna' if OPTIONS.DNA else 'amino')

    alignment_basename = '%s_%s_%s' % (ab_basename, splitext(basename(OPTIONS.NEUT_SQLITE3_DB))[0], __VERSION__)

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords = collect_seqrecords_from_db(OPTIONS.NEUT_SQLITE3_DB, antibody, dna=OPTIONS.DNA)
    alignment, refseq_offs = generate_alignment(seqrecords, alignment_basename, is_HXB2, OPTIONS)
    colfilter = None
    if OPTIONS.PHYLOFILTER:
        colfilter = PhyloFilter(
            alph,
            _PHYLOFILTER_BATCHFILE,
            is_HXB2,
            lambda x: False
        )
        colnames, x = colfilter.learn(alignment, refseq_offs)
    else:
        colfilter = NaiveFilter(
            alph,
            OPTIONS.MAX_CONSERVATION,
            OPTIONS.MIN_CONSERVATION,
            OPTIONS.MAX_GAP_RATIO,
            is_HXB2,
            lambda x: False # TODO: add the appropriate filter function based on the args here
        )

        if sim is None:
            colnames, x = colfilter.learn(alignment, refseq_offs)
            # TODO: I don't think we need to binarize the colnames here, though we can if we want.
            # I need to think more about how to properly handle this case.
        else:
            if OPTIONS.SIM_EPI_N is None:
                OPTIONS.SIM_EPI_N = len(seqrecords)

    # compute features
    for target in OPTIONS.TARGETS:

        forward_initval = 1 if OPTIONS.FORWARD_SELECT else OPTIONS.NUM_FEATURES
        forward_select = None
        results = None
        for num_features in range(forward_initval, OPTIONS.NUM_FEATURES + 1):

            if sim is None:
                yextractor = ClassExtractor(
                    seqrecord_to_ic50s,
                    lambda row: is_HXB2(row) or False, # TODO: again filtration function
                    lambda x: x < OPTIONS.IC50LT if target == 'lt' else x > OPTIONS.IC50GT
                )
                y = yextractor.extract(alignment)

#             if OPTIONS.PHYLOFILTER:
#                 np.savez('phylo.y.npz', {'data': y})
#             else:
#                 with open(antibody + '.mrmr', 'w') as fh:
#                     print('class,' + ','.join(colnames), file=fh)
#                     print('\n'.join(','.join('%d' % v for v in ([y[i]] + (x[i, :].tolist()))) for i in range(x.shape[0])), file=fh)

            # simulations, ho!
            for i in range(sim.runs if sim is not None else 1):

                # here is where the sequences must be generated for the random sequence and random epitope simulations
                if sim is not None:
                    alignment = sim.generate_sequences(
                        N=OPTIONS.SIM_EPI_N,
                        idfmt='%s|||',
                        noise=OPTIONS.SIM_EPI_NOISE,
                        mutation_rate=OPTIONS.SIM_EPI_MUT_RATE,
                        alphabet=alph
                    )

                    colnames, x = colfilter.learn(alignment, {}) # XXX: refseq_offs needed here?

                    # simulates the epitope and assigns the appropriate class
                    epi_def = sim.simulate_epitope(
                        alignment,
                        alph,
                        colnames,
                        OPTIONS.SIM_EPI_SIZE,
                        OPTIONS.SIM_EPI_PERCENTILE,
                    )

                    if epi_def is not None:
                        print('********************* SIMULATED EPITOPE DESCRIPTION (%d) *********************\n' % OPTIONS.SIM_EPI_SIZE, file=sys.stdout)
                        print('%s\n' % str(epi_def), file=sys.stdout)

                optstat = DiscretePerfStats.MINSTAT
                if OPTIONS.ACCURACY:
                    optstat = DiscretePerfStats.ACCURACY
                elif OPTIONS.PPV:
                    optstat = DiscretePerfStats.PPV
                elif OPTIONS.NPV:
                    optstat = DiscretePerfStats.NPV
                elif OPTIONS.SENSITIVITY:
                    optstat = DiscretePerfStats.SENSITIVITY
                elif OPTIONS.SPECIFICITY:
                    optstat = DiscretePerfStats.SPECIFICITY
                elif OPTIONS.FSCORE:
                    optstat = DiscretePerfStats.FSCORE

                C_begin, C_end, C_step = OPTIONS.LOG2C
                recip = 1
                if isinstance(C_step, float):
                    recip = 1. / C_step
                    C_begin, C_end = int(recip * C_begin), int(recip * C_end)
                    C_step = 1
                C_range = [pow(2., float(C) / recip) for C in range(C_begin, C_end + 1, C_step)]

                if OPTIONS.MRMR_NORMALIZE:
                    DiscreteMrmr._NORMALIZED = True

                crossvalidator = SelectingNestedCrossValidator(
                    classifier_cls=LinearSvm,
                    selector_cls=PhyloMrmr if OPTIONS.PHYLOFILTER else DiscreteMrmr,
                    folds=OPTIONS.CV_FOLDS,
                    gridsearch_kwargs={ 'C': C_range },
                    classifier_kwargs={},
                    selector_kwargs={
                        'num_features': num_features,
                        'method': DiscreteMrmr.MAXREL if OPTIONS.MAXREL else \
                                  DiscreteMrmr.MID if OPTIONS.MRMR_METHOD == 'MID' else \
                                  DiscreteMrmr.MIQ
                    },
                    validator_cls=CrossValidator,
                    validator_kwargs={
                        'folds': OPTIONS.CV_FOLDS-1,
                        'scorer_cls': DiscretePerfStats,
                        'scorer_kwargs': { 'optstat': optstat }
                    },
                    scorer_cls=DiscretePerfStats,
                    scorer_kwargs={ 'optstat': optstat },
                    weights_func='weights' # we MUST specify this or it will be set to lambda: None
                )

                new_results = crossvalidator.crossvalidate(x, y, classifier_kwargs={}, extra=extract_feature_weights)

                if results is not None and new_results.stats.get() <= results.stats.get():
                    break

                results = new_results
                forward_select = num_features

        # the alignment reflects the number of sequences either naturally,
        # or through SIM_EPI_N, which reflects the natural number anyway.
        meta = make_output_meta(OPTIONS, len(alignment), target, antibody, forward_select)
        ret = cv_results_to_output(results, colnames, meta)

        if isinstance(OPTIONS.OUTPUT, str):
            with open(OPTIONS.OUTPUT, 'w') as fh:
                print(pretty_fmt_results(ret), file=fh)
        else:
            print(pretty_fmt_results(ret))

    return 0


if __name__ == '__main__':
    sys.exit(main())
