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

import sqlite3, sys

from codecs import getwriter
from math import ceil, copysign, log10, sqrt
from operator import itemgetter
from optparse import OptionParser
from os import close, remove, rename
from os.path import basename, dirname, exists, join, realpath, splitext
from random import gauss, random, seed
from re import sub, match
from tempfile import mkstemp

# add biopython, cvxopt, pil to my path if I'm on Mac OS X
import platform
if platform.system().lower() == 'darwin':
    __subpath = join('python%d%d' % sys.version_info[:2], 'lib', 'python')
    # numpy must go first
    for module in ('numpy-1.5.1', 'biopython-1.56', 'cvxopt-1.1.3', 'cython-0.14.1', 'mlpy-2.2.2', 'pil-1.1.7'):
        path = join(dirname(realpath(__file__)), 'contrib', module, __subpath)
        if exists(path):
            sys.path.insert(0, path)

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from idepi import *

from numpy import mean, median, std

__version__ = 0.5

_WORKING_DIR      = dirname(realpath(__file__)) # '/Users/Lance/Pond'
_HXB2_DNA_FASTA   = join(_WORKING_DIR, 'res', 'hxb2_dna.fa')
_HXB2_AMINO_FASTA = join(_WORKING_DIR, 'res', 'hxb2_pep.fa')

_DEFAULT_NUM_FEATURES = 10

_RAND_SEQ = 'randseq'
_RAND_TARGET = 'randtarget'
_RAND_EPI = 'randepi'
_RAND_DUMB = 'randdumbepi'
_SIM_VALS = (_RAND_SEQ, _RAND_TARGET, _RAND_EPI, _RAND_DUMB)
_RAND_SEQ_STOCKHOLM = join(_WORKING_DIR, 'res', 'randhivenvpep_final.sto')

_PHYLOFILTER_BATCHFILE = join(_WORKING_DIR, 'res', 'hyphy', 'CorrectForPhylogeny.bf')

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

_TEST_AMINO_MRMR = ['1,1,0,1,0,1,1,0,1,0,1,0,1,0,1',
                    '0,1,0,1,0,1,1,0,1,0,1,0,1,1,0',
                    '0,1,0,1,0,1,0,1,0,1,1,0,1,1,0',
                    '1,0,1,0,1,1,1,0,1,0,1,1,0,0,1']

_TEST_STANFEL_MRMR = ['1,1,0,1,1,1,0,1,0,1,0,1,0,1',
                      '0,1,0,1,1,1,0,1,0,1,0,1,1,0',
                      '0,1,0,1,1,0,1,0,1,1,0,1,1,0',
                      '1,0,1,1,1,1,0,1,0,1,1,0,0,1']

_TEST_AMINO_SVM = ['1 1:1 3:1 5:1 6:1 8:1 10:1 12:1 14:1',
                   '0 1:1 3:1 5:1 6:1 8:1 10:1 12:1 13:1',
                   '0 1:1 3:1 5:1 7:1 9:1 10:1 12:1 13:1',
                   '1 2:1 4:1 5:1 6:1 8:1 10:1 11:1 14:1']

_TEST_STANFEL_SVM = ['1 1:1 3:1 4:1 5:1 7:1 9:1 11:1 13:1',
                     '0 1:1 3:1 4:1 5:1 7:1 9:1 11:1 12:1',
                     '0 1:1 3:1 4:1 6:1 8:1 9:1 11:1 12:1',
                     '1 2:1 3:1 4:1 5:1 7:1 9:1 10:1 13:1']


OPTIONS = None


def optparse_extend(option, opt_str, value, parser):
    if getattr(parser.values, option.dest, None) is None:
        setattr(parser.values, option.dest, list())
    getattr(parser.values, option.dest).extend(value)


def optparse_csv(option, opt_str, value, parser):
    setattr(parser.values, option.dest, value.split(','))


def setup_option_parser():

    parser = OptionParser(usage = '%prog [options] ANTIBODY')

    #                 option            action = 'store'        callback                    type             nargs = 1  dest
    parser.add_option('--hmmalign',                                                         type = 'string',            dest = 'HMMER_ALIGN_BIN')
    parser.add_option('--hmmbuild',                                                         type = 'string',            dest = 'HMMER_BUILD_BIN')
    parser.add_option('--hmmiter',                                                          type = 'int',               dest = 'HMMER_ITER')
    parser.add_option('--mrmr',         action = 'store_false',                                                         dest = 'MAXREL')
    parser.add_option('--mrmrprog',                                                         type = 'string',            dest = 'MRMR_BIN')
    parser.add_option('--mrmrmethod',                                                       type = 'string',            dest = 'MRMR_METHOD')
    parser.add_option('--maxrel',       action = 'store_true',                                                          dest = 'MAXREL')
    parser.add_option('--filter',       action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'FILTER')
    parser.add_option('--numfeats',                                                         type = 'int',               dest = 'NUM_FEATURES')
    parser.add_option('--subtypes',     action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'SUBTYPES')
    parser.add_option('--svmtrain',                                                         type = 'string',            dest = 'SVM_TRAIN_BIN')
    parser.add_option('--svmpredict',                                                       type = 'string',            dest = 'SVM_PREDICT_BIN')
    parser.add_option('--svmargs',      action = 'callback',    callback = optparse_extend, type = 'string', nargs = 2, dest = 'SVM_ARGS')
    parser.add_option('--nested',       action = 'store_true',                                                          dest = 'NESTED')
    parser.add_option('--log2c',        action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'LOG2C')
    parser.add_option('--weighting',    action = 'store_true',                                                          dest = 'WEIGHTING')
    parser.add_option('--accuracy',     action = 'store_true',                                                          dest = 'ACCURACY')
    parser.add_option('--ppv',          action = 'store_true',                                                          dest = 'PPV')
    parser.add_option('--precision',    action = 'store_true',                                                          dest = 'PPV')
    parser.add_option('--npv'     ,     action = 'store_true',                                                          dest = 'NPV')
    parser.add_option('--sensitivity',  action = 'store_true',                                                          dest = 'SENSITIVITY')
    parser.add_option('--recall',       action = 'store_true',                                                          dest = 'SENSITIVITY')
    parser.add_option('--specificity',  action = 'store_true',                                                          dest = 'SPECIFICITY')
    parser.add_option('--tnr',          action = 'store_true',                                                          dest = 'SPECIFICITY')
    parser.add_option('--fscore',       action = 'store_true',                                                          dest = 'FSCORE')
    parser.add_option('--amino',        action = 'store_true',                                                          dest = 'AMINO')
    parser.add_option('--dna',          action = 'store_true',                                                          dest = 'DNA')
    parser.add_option('--stanfel',      action = 'store_true',                                                          dest = 'STANFEL')
    parser.add_option('--cv',                                                               type = 'int',               dest = 'CV_FOLDS')
    parser.add_option('--loocv',        action = 'store_true',                                                          dest = 'LOOCV')
    parser.add_option('--targets',      action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'TARGETS')
    parser.add_option('--maxcon',                                                           type = 'float',             dest = 'MAX_CONSERVATION')
    parser.add_option('--maxgap',                                                           type = 'float',             dest = 'MAX_GAP_RATIO')
    parser.add_option('--mincon',                                                           type = 'float',             dest = 'MIN_CONSERVATION')
    parser.add_option('--ic50gt',                                                           type = 'float',             dest = 'IC50GT')
    parser.add_option('--ic50lt',                                                           type = 'float',             dest = 'IC50LT')
    parser.add_option('--neuts',                                                            type = 'string',            dest = 'NEUT_SQLITE3_DB')
    parser.add_option('--hxb2',                                                             type = 'string',            dest = 'HXB2_FASTA')
    parser.add_option('--ids',          action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'HXB2_IDS')
    parser.add_option('--test',         action = 'store_true',                                                          dest = 'TEST')
    parser.add_option('--sim',                                                              type = 'string',            dest = 'SIM')
    parser.add_option('--simruns',                                                          type = 'int',               dest = 'SIM_RUNS')
    parser.add_option('--simepisize',                                                       type = 'int',               dest = 'SIM_EPI_SIZE')
    parser.add_option('--simepimutrate',                                                    type = 'float',             dest = 'SIM_EPI_MUT_RATE')
    parser.add_option('--simepiseqnum',                                                     type = 'int',               dest = 'SIM_EPI_N')
    parser.add_option('--simepinoise',                                                      type = 'float',             dest = 'SIM_EPI_NOISE')
    parser.add_option('--simepiperc',                                                       type = 'float',             dest = 'SIM_EPI_PERCENTILE')
    parser.add_option('--seed',                                                             type = 'int',               dest = 'RAND_SEED')
    parser.add_option('--phylofilt',    action = 'store_true',                                                          dest = 'PHYLOFILTER')

    parser.set_defaults(HMMER_ALIGN_BIN    = join(_WORKING_DIR, 'contrib', 'hmmer-3.0', 'src', 'hmmalign'))
    parser.set_defaults(HMMER_BUILD_BIN    = join(_WORKING_DIR, 'contrib', 'hmmer-3.0', 'src', 'hmmbuild'))
    parser.set_defaults(HMMER_ITER         = 8)
    parser.set_defaults(MRMR_BIN           = join(_WORKING_DIR, 'contrib', 'mrmr_c_src', 'mrmr'))
    parser.set_defaults(MRMR_METHOD        = 'MID')
    parser.set_defaults(MAXREL             = False)
    parser.set_defaults(FILTER             = [])
    parser.set_defaults(NUM_FEATURES       = -1)
    parser.set_defaults(SUBTYPES           = [])
    parser.set_defaults(SVM_TRAIN_BIN      = join(_WORKING_DIR, 'contrib', 'libsvm-3.0', 'svm-train'))
    parser.set_defaults(SVM_PREDICT_BIN    = join(_WORKING_DIR, 'contrib', 'libsvm-3.0', 'svm-predict'))
    parser.set_defaults(SVM_ARGS           = ['-t', '0', '-h', '0'])
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
    parser.set_defaults(TARGETS            = ['gt', 'lt'])
    parser.set_defaults(MAX_CONSERVATION   = 1. ) # 93.)
    parser.set_defaults(MAX_GAP_RATIO      = 0.1 ) # 93.)
    parser.set_defaults(MIN_CONSERVATION   = 1. ) # 33.)
    parser.set_defaults(IC50GT             = 20.)
    parser.set_defaults(IC50LT             = 2.)
    parser.set_defaults(NEUT_SQLITE3_DB    = join(_WORKING_DIR, 'res', 'allneuts.sqlite3'))
    parser.set_defaults(HXB2_FASTA         = _HXB2_AMINO_FASTA)
    parser.set_defaults(HXB2_IDS           = ['9629357', '9629363'])
    parser.set_defaults(SIM                = '') # can be 'randtarget' for now
    parser.set_defaults(SIM_RUNS           = 1) # can be 'randtarget' for now
    parser.set_defaults(SIM_EPI_SIZE       = 10)
    parser.set_defaults(SIM_EPI_MUT_RATE   = 0.01)
    parser.set_defaults(SIM_EPI_N          = None) # default is to use len(abrecords)
    parser.set_defaults(SIM_EPI_NOISE      = 0.08)
    parser.set_defaults(SIM_EPI_PERCENTILE = 0.5)
    parser.set_defaults(RAND_SEED          = 42) # magic number for determinism
    parser.set_defaults(PHYLOFILTER        = False)

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
        print >> fh, _TEST_AMINO_STO
        fh.close()

        seq_table = SeqTable(sto_filename, SeqTable.AMINO_ALPHABET)

        for OPTIONS.STANFEL in (True, False):

            if OPTIONS.STANFEL:
                OPTIONS.AMINO = False
                _TEST_NAMES = _TEST_STANFEL_NAMES
                _TEST_MRMR = _TEST_STANFEL_MRMR
                _TEST_SVM = _TEST_STANFEL_SVM
            else:
                OPTIONS.AMINO = True
                _TEST_NAMES = _TEST_AMINO_NAMES
                _TEST_MRMR = _TEST_AMINO_MRMR
                _TEST_SVM = _TEST_AMINO_SVM

            alph = Alphabet(Alphabet.STANFEL if OPTIONS.STANFEL else Alphabet.DNA if OPTIONS.DNA else Alphabet.AMINO)

            all_feature_names = generate_feature_names(seq_table, alph)

            # TODO: fix this stupidity
            feature_idxs, feature_names = compute_relevant_features_and_names(
                seq_table,
                alph,
                all_feature_names,
                {
                    'max': OPTIONS.MAX_CONSERVATION,
                    'min': OPTIONS.MIN_CONSERVATION,
                    'gap': OPTIONS.MAX_GAP_RATIO
                },
                OPTIONS.FILTER,
            )

            feature_names_ = list(feature_names)

            # test the feature names portion
            try:
                assert(len(feature_names_) == len(_TEST_NAMES))
            except AssertionError, e:
                print 'gen:   %s\ntruth: %s' % (feature_names_, _TEST_NAMES)
                raise e

            mrmr_header = ','.join(['class'] + _TEST_NAMES)

            for name in _TEST_NAMES:
                try:
                    assert(name in feature_names_)
                    del feature_names_[feature_names_.index(name)]
                except AssertionError, e:
                    print >> sys.stderr, 'ERROR: \'%s\' not found in %s' % (name, ', '.join(feature_names_))
                    raise e

            # test mRMR and LSVM file generation
            for target in OPTIONS.TARGETS:
                data = generate_relevant_data(
                    feature_names,
                    seq_table,
                    Alphabet(mode=Alphabet.AMINO),
                    feature_idxs,
                    target
                )

                x, y = data.tondarrays()

                # generate and test the mRMR portion
                mrmr = Mrmr(num_features=OPTIONS.NUM_FEATURES, method=Mrmr.MAXREL if OPTIONS.MAXREL else Mrmr.MID if OPTIONS.MRMR_METHOD == 'MID' else Mrmr.MIQ)

                mrmr.select(x, y)

                x = mrmr.subset(x)

                lsvm = LinearSvm()
                lsvm.learn(x, y)

    finally:
        remove(sto_filename)

    print >> sys.stderr, 'ALL TESTS PASS'


def fix_hxb2_fasta():
    '''If DNA mode was selected but the AMINO reference sequence is still in place, fix it'''
    if OPTIONS.DNA == True and OPTIONS.HXB2_FASTA == _HXB2_AMINO_FASTA:
        OPTIONS.HXB2_FASTA = _HXB2_DNA_FASTA


def generate_alignment(seqrecords, filename):
    if not exists(filename) and OPTIONS.SIM != _RAND_DUMB:
        generate_alignment_from_SeqRecords(
            OPTIONS.HXB2_FASTA,
            seqrecords,
            OPTIONS.HMMER_ALIGN_BIN,
            OPTIONS.HMMER_BUILD_BIN,
            OPTIONS.HMMER_ITER,
            filename,
            dna=True if OPTIONS.DNA else False
        )
    elif OPTIONS.SIM == _RAND_DUMB:
        # we're assuming pre-aligned because they're all generated from the same refseq
        fh = open(filename, 'w')
        SeqIO.write(seqrecords, fh, 'stockholm')
        fh.close()

    with open(filename) as fh:
        alignment = AlignIO.read(fh, 'stockholm')

    return alignment


def collect_SeqTable_and_feature_names(alignment_filename, seq_records):

    skip_func = None
    if len(OPTIONS.SUBTYPES) != 0:
        skip_func = lambda x: len([i for i in x.split('|')[1] if i in tuple(OPTIONS.SUBTYPES)]) <= 0

    alph = Alphabet(Alphabet.STANFEL if OPTIONS.STANFEL else Alphabet.DNA if OPTIONS.DNA else Alphabet.AMINO)

    seq_table = SeqTable(alignment_filename, SeqTable.DNA_ALPHABET if OPTIONS.DNA else SeqTable.AMINO_ALPHABET, is_HXB2, skip_func)

    all_feature_names = generate_feature_names(seq_table, alph)

    # set the right number of CV_FOLDS for Leave-One-Out-Crossvalidation
    if OPTIONS.LOOCV:
        OPTIONS.CV_FOLDS = len([r for r in seq_table.rows() if not is_HXB2(r.id)])

    # parition the table into CV_FOLDS groups for cross-validation
    seq_table.partition(OPTIONS.CV_FOLDS)

    return seq_table, all_feature_names


def main(argv = sys.argv):
    global OPTIONS

    # fix some stupid bugs
    sys.stdout = getwriter('utf8')(sys.stdout)

    # so some option parsing
    option_parser = setup_option_parser()
    (OPTIONS, args) = option_parser.parse_args(argv)

    # do some argument parsing
    if OPTIONS.TEST:
        run_tests()
        return 0

    if OPTIONS.SIM != '' and OPTIONS.SIM not in _SIM_VALS:
        option_parser.error('option --sim takes one of %s' % ', '.join(_SIM_VALS))

    if OPTIONS.RAND_SEED is not None:
        seed(OPTIONS.RAND_SEED)

    if len(args) != 2:
        option_parser.error('ANTIBODY is a required argument')

    if not set(OPTIONS.TARGETS).issubset(set(['lt', 'gt'])):
        option_parser.error('option --targets takes either or both: lt gt')

    if not OPTIONS.MRMR_METHOD in ('MIQ', 'MID'):
        option_parser.error('option --mrmrmethod takes either MIQ or MID')

    # check to make sure our mode is exclusive, and set the default (AMINO) if none is set
    if sum([1 for v in (OPTIONS.AMINO, OPTIONS.DNA, OPTIONS.STANFEL) if v]) > 1:
        option_parser.error('options --amino, --dna, and --stanfel are mutually exclusive')
    elif sum([1 for v in (OPTIONS.AMINO, OPTIONS.DNA, OPTIONS.STANFEL) if v]) == 0:
        OPTIONS.AMINO = True

    if sum([1 for v in (OPTIONS.ACCURACY, OPTIONS.PPV, OPTIONS.NPV, OPTIONS.SENSITIVITY, OPTIONS.SPECIFICITY, OPTIONS.FSCORE) if v]) > 1:
        option_parser.error('options --accuracy, --ppv/--precision, --npv, --sensitivity/--recall, --specificity/--tnr, --fscore are mutually exclusive')

    try:
        if len(OPTIONS.LOG2C) != 3:
            raise ValueError
        OPTIONS.LOG2C = [int(OPTIONS.LOG2C[0]), int(OPTIONS.LOG2C[1]), float(OPTIONS.LOG2C[2])]
    except ValueError, e:
        option_parser.error('option --log2c takes an argument of the form C_BEGIN,C_END,C_STEP')

    # validate the antibody argument, currently a hack exists to make PG9/PG16 work
    # TODO: Fix pg9/16 hax
    antibody = args[1]
    valid_antibodies = sorted(get_valid_antibodies_from_db(OPTIONS.NEUT_SQLITE3_DB), key = lambda x: x.strip())
    if antibody not in valid_antibodies:
        if ' ' + antibody not in valid_antibodies:
            option_parser.error('%s not in the list of permitted antibodies: \n  %s' % (antibody, '\n  '.join([ab.strip() for ab in valid_antibodies])))
        else:
            antibody = ' ' + antibody

    # validate the subtype option
    valid_subtypes = sorted(get_valid_subtypes_from_db(OPTIONS.NEUT_SQLITE3_DB), key = lambda x: x.strip().upper())
    for subtype in OPTIONS.SUBTYPES:
        if subtype not in valid_subtypes:
            option_parser.error('%s not in the list of permitted subtypes: \n  %s' % (subtype, '\n  '.join([st.strip() for st in valid_subtypes])))

    if OPTIONS.SIM in (_RAND_EPI, _RAND_SEQ) and OPTIONS.DNA:
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
            print >> sys.stderr, 'warning: clamping --numfeats to sizeof(--filter) = %d' % OPTIONS.NUM_FEATURES
    else: # len(OPTIONS.FILTER) == 0
        if OPTIONS.NUM_FEATURES == -1:
            OPTIONS.NUM_FEATURES = _DEFAULT_NUM_FEATURES

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
        if OPTIONS.SIM == _RAND_DUMB:
            hxb2fh = open(OPTIONS.HXB2_FASTA)
            for record in SeqIO.parse(hxb2fh, 'fasta'):
                hxb2seq = str(record.seq)
                break
            hxb2fh.close()
            sim = DumbSimulation(OPTIONS.SIM_RUNS, Simulation.EPITOPE, hxb2seq)
        elif OPTIONS.SIM == _RAND_EPI or OPTIONS.SIM == _RAND_SEQUENCE:
            sim = MarkovSimulation(OPTIONS.SIM_RUNS, Simulation.SEQUENCE if OPTIONS.SIM != _RAND_EPI else Simulation.EPITOPE, _RAND_SEQ_STOCKHOLM)
        else:
            raise ValueError('Unknown simulation type `%s\'' % OPTIONS.SIM)

    ab_basename = '%s%s_%s' % (antibody, '_randseq' if sim is not None and sim.mode == Simulation.SEQUENCE else '', 'dna' if OPTIONS.DNA else 'amino')

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    abrecords = collect_AbRecords_from_db(OPTIONS.NEUT_SQLITE3_DB, antibody)

    # do not generate the sequences here for the random sequence or random epitope simulations
    alignment_filename = None
    if sim is None:
        alignment_filename = '%s_%s_%s.sto' % (ab_basename, splitext(basename(OPTIONS.NEUT_SQLITE3_DB))[0], __version__)

        # generate an alignment using HMMER if it doesn't already exist
        seqrecords = [r.to_SeqRecord(dna=True if OPTIONS.DNA else False) for r in abrecords]
        alignment = generate_alignment(seqrecords, alignment_filename)
        colfilter = None
        if OPTIONS.PHYLOFILTER:
            colfilter = PhyloFilter(
                alph,
                _PHYLOFILTER_BATCHFILE,
                is_HXB2,
                lambda x: False
            )
        else:
            colfilter = NaiveFilter(
                alph,
                OPTIONS.MAX_CONSERVATION,
                OPTIONS.MIN_CONSERVATION,
                OPTIONS.MAX_GAP_RATIO,
                is_HXB2,
                lambda x: False # TODO: add the appropriate filter function based on the args here
            )
        colnames, x = colfilter.filter(alignment)
    else:
        if OPTIONS.SIM_EPI_N is None:
            OPTIONS.SIM_EPI_N = len(abrecords)


    # compute features
    for target in OPTIONS.TARGETS:
        results = None

        yextractor = ClassExtractor(
            id_to_float,
            lambda row: is_HXB2(row) or False, # TODO: again filtration function
            lambda x: x < OPTIONS.IC50LT if target == 'lt' else x > OPTIONS.IC50GT
        )

        y = yextractor.extract(alignment)

        # simulations, ho!
        for i in xrange(sim.runs if sim is not None else 1):

            # here is where the sequences must be generated for the random sequence and random epitope simulations
            if sim is not None:
                alignment_filename = '%s_%s_%d_%s.sto' % (ab_basename, OPTIONS.SIM, i + 1, __version__)
                seq_records = sim.generate_sequences(
                    alignment_filename,
                    N=OPTIONS.SIM_EPI_N,
                    idfmt='%s|||',
                    noise=OPTIONS.SIM_EPI_NOISE,
                    mutation_rate=OPTIONS.SIM_EPI_MUT_RATE,
                    alphabet=alph
                )

                seq_table, all_feature_names = collect_SeqTable_and_feature_names(alignment_filename, seq_records)

                # simulates the epitope and assigns the appropriate class
                epi_def = sim.simulate_epitope(
                    seq_table,
                    alph,
                    all_feature_names,
                    OPTIONS.SIM_EPI_SIZE,
                    OPTIONS.SIM_EPI_PERCENTILE,
                )

                if epi_def is not None:
                    print >> sys.stdout, '********************* SIMULATED EPITOPE DESCRIPTION (%d) *********************\n' % OPTIONS.SIM_EPI_SIZE
                    print >> sys.stdout, '%s\n' % str(epi_def)

            optstat = PerfStats.MINSTAT
            if OPTIONS.ACCURACY:
                optstat = PerfStats.ACCURACY
            elif OPTIONS.PPV:
                optstat = PerfStats.PPV
            elif OPTIONS.NPV:
                optstat = PerfStats.NPV
            elif OPTIONS.SENSITIVITY:
                optstat = PerfStats.SENSITIVITY
            elif OPTIONS.SPECIFICITY:
                optstat = PerfStats.SPECIFICITY
            elif OPTIONS.FSCORE:
                optstat = PerfStats.FSCORE

            C_begin, C_end, C_step = OPTIONS.LOG2C
            recip = 1
            if isinstance(C_step, float):
                recip = 1. / C_step
                C_begin, C_end = int(recip * C_begin), int(recip * C_end)
                C_step = 1
            C_range = [pow(2., float(C) / recip) for C in xrange(C_begin, C_end + 1, C_step)]

            crossvalidator = SelectingNestedCrossValidator(
                classifiercls=LinearSvm,
                selectorcls=FastCaimMrmr if OPTIONS.PHYLOFILTER else DiscreteMrmr,
                folds=OPTIONS.CV_FOLDS,
                cv={},
                optstat=optstat,
                gs={ 'C': C_range },
                fs={
                    'num_features': OPTIONS.NUM_FEATURES,
                    'method': DiscreteMrmr.MAXREL if OPTIONS.MAXREL else \
                              DiscreteMrmr.MID if OPTIONS.MRMR_METHOD == 'MID' else \
                              DiscreteMrmr.MIQ
                    }
            )

            results = crossvalidator.crossvalidate(x, y, cv={}, extra=lambda x: { 'features': x.features(), 'weights': x.classifier.weights() })

            # remove the alignment
            if OPTIONS.SIM in (_RAND_DUMB, _RAND_EPI):
                remove(alignment_filename)

        # print stats and results:
#         print >> sys.stdout, '********************* REPORT FOR ANTIBODY %s IC50 %s *********************' % \
#           (antibody, '< %d' % OPTIONS.IC50LT if target == 'lt' else '> %d' % OPTIONS.IC50GT)
#
#         if OPTIONS.SIM not in (_RAND_SEQ, _RAND_TARGET, _RAND_EPI):
#             fmt = ('', OPTIONS.CV_FOLDS, '')
#         else:
#             fmt = ('%d-run ' % OPTIONS.SIM_RUNS, OPTIONS.CV_FOLDS, ' per run')

        statsdict = results['stats'].todict()

        featureweights = {}
        for i in xrange(len(results['extra'])):
            assert(len(results['extra'][i]['features']) >= len(results['extra'][i]['weights']))
            for j in xrange(len(results['extra'][i]['features'])):
                v = results['extra'][i]['weights'][j] if j < len(results['extra'][i]['weights']) else 0.
                k = results['extra'][i]['features'][j]
                if k not in featureweights:
                    featureweights[k] = []
                featureweights[k].append(int(copysign(1, v)))

        weightsdict = {}
        for idx, weights in featureweights.items():
            val = NormalValue(int, weights)
            weightsdict[colnames[idx]] = val

        # remove minstat 'cause we don't want it here..
        if 'Minstat' in statsdict:
            del statsdict['Minstat']

        ret = {}

        ret['statistics'] = dict([(k.lower(), { 'mean': v.mu, 'std': sqrt(v.sigma) }) for k, v in statsdict.items()])
        ret['weights'] = [{ 'position': k, 'value': { 'mean': v.mu, 'std': sqrt(v.sigma), 'N': len(v) } } for k, v in sorted(
            weightsdict.items(),
            key=lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x[0]))
        )]

        print >> sys.stdout, '{\n  "statistics": {'#  % OPTIONS.CV_FOLDS

        stat_len = max([len(k) for k in ret['statistics'].keys()]) + 3
        mean_len = max([len('%.6f' % v['mean']) for v in ret['statistics'].values()])
        std_len = max([len('%.6f' % v['std']) for v in ret['statistics'].values()])
        fmt = u'{ "mean": %%%d.6f, "std": %%%d.6f }' % (mean_len, std_len)
        output = [u'    %-*s %s' % (
            stat_len, '"%s":' % k,
            fmt % (v['mean'], v['std'])
        ) for k, v in sorted(ret['statistics'].items(), key=itemgetter(0))]
        print >> sys.stdout, ',\n'.join(output)

        print >> sys.stdout, '  },\n  "weights": ['

        if len(ret['weights']) > 0:
            name_len = max([len(v['position']) for v in ret['weights']]) + 3
            mean_len = max([len('% .6f' % v['value']['mean']) for v in ret['weights']])
            std_len = max([len('%.6f' % v['value']['std']) for v in ret['weights']])
            N_len = max([len('%d' % v['value']['N']) for v in ret['weights']])
            fmt = u'{ "mean": %%%d.6f, "std": %%%d.6f, "N": %%%dd }' % (mean_len, std_len, N_len)
            output = [u'    { "position": %-*s "value": %s }' % (
                name_len, u'"%s",' % v['position'],
                fmt % (
                    v['value']['mean'],
                    v['value']['std'],
                    v['value']['N']
                ),
            ) for v in sorted(ret['weights'], key=lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x['position'])))]
            print >> sys.stdout, ',\n'.join(output)

        print >> sys.stdout, '  ]\n}'

    return 0


if __name__ == '__main__':
    sys.exit(main())
