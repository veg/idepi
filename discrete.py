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
from os import remove, rename
from os.path import basename, dirname, exists, join, realpath, splitext
from random import gauss, randint, random, seed
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

_STANFEL = {'A': 0,
            'C': 0,
            'G': 0,
            'I': 0,
            'L': 0,
            'M': 0,
            'P': 0,
            'S': 0,
            'T': 0,
            'V': 0,
            'D': 1,
            'E': 1,
            'N': 1,
            'Q': 1,
            'F': 2,
            'W': 2,
            'Y': 2,
            'H': 3,
            'K': 3,
            'R': 3,
            'X': 4,
            '-': 5} # this is a space, but we cant have 'keyword be an expression'

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
    parser.set_defaults(SIM_EPI_N          = None) # default is to use len(ab_records)
    parser.set_defaults(SIM_EPI_NOISE      = 0.08)
    parser.set_defaults(SIM_EPI_PERCENTILE = 0.5)
    parser.set_defaults(RAND_SEED          = 42) # magic number for determinism

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

    sto_filename = mkstemp()[1]

    try:
        fh = open(sto_filename, 'w')
        print >> fh, _TEST_AMINO_STO
        fh.close()

        seq_table = SeqTable(sto_filename, AMINO_ALPHABET)

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

            all_feature_names = generate_feature_names(seq_table)

            # TODO: fix this stupidity
            feature_idxs, feature_names = compute_relevant_features_and_names(seq_table, all_feature_names)

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
                data = generate_relevant_data(seq_table, feature_names, feature_idxs, target)

                x, y = data.tondarrays()

                # generate and test the mRMR portion
                mrmr = Mrmr(num_features=OPTIONS.NUM_FEATURES, method=MAXREL if OPTIONS.MAXREL else MID if OPTIONS.MRMR_METHOD == 'MID' else MIQ)

                mrmr.select(x, y)

                x = mrmr.subset(x)

                lsvm = LinearSvm()
                lsvm.learn(x, y)

    finally:
        remove(sto_filename)

    print >> sys.stderr, 'ALL TESTS PASS'


def get_valid_subtypes_from_db():
    conn = sqlite3.connect(OPTIONS.NEUT_SQLITE3_DB)
    curr = conn.cursor()

    curr.execute('''select distinct SUBTYPE from GENO_REPORT''')

    valid_subtypes = [r[0] for r in curr if r[0].strip() != '']

    conn.close()

    return valid_subtypes


def get_valid_antibodies_from_db():
    conn = sqlite3.connect(OPTIONS.NEUT_SQLITE3_DB)
    curr = conn.cursor()

    curr.execute('''select distinct ANTIBODY from NEUT''')

    valid_antibodies = [r[0] for r in curr]

    conn.close()

    return valid_antibodies


def fix_hxb2_fasta():
    '''If DNA mode was selected but the AMINO reference sequence is still in place, fix it'''
    if OPTIONS.DNA == True and OPTIONS.HXB2_FASTA == _HXB2_AMINO_FASTA:
        OPTIONS.HXB2_FASTA = _HXB2_DNA_FASTA


class ABRecord(object):

    def __init__(self, id, seq, subtype, ab, ic50):
        self.id = id
        self.dna_seq, self.amino_seq = OrfList(seq)[0]
        self.subtype = subtype.upper()
        self.antibody = ab
        ic50 = ic50.strip()
        if '>' in ic50:
            ic50 = 25.
        else:
            ic50 = float(ic50)
        self.ic50 = ic50

    def to_SeqRecord(self):
        if OPTIONS.DNA:
            _seq = self.dna_seq
        else:
            _seq = self.amino_seq
        return SeqRecord(_seq, '%s|%s|%s|%s' % (self.id, self.subtype, self.antibody, self.ic50))


def collect_ABRecords_from_db(antibody):
    conn = sqlite3.connect(OPTIONS.NEUT_SQLITE3_DB)
    curr = conn.cursor()

    curr.execute('''
    select distinct S.ID as ID, S.SEQ as SEQ, G.SUBTYPE as SUBTYPE, N.AB as AB, N.IC50 as IC50 from
    (select ACCESSION_ID as ID, RAW_SEQ as SEQ from SEQUENCE group by ACCESSION_ID) as S join
    (select ACCESSION_ID as ID, SUBTYPE from GENO_REPORT group by ACCESSION_ID) as G join
    (select ACCESSION_ID as ID, ANTIBODY as AB, IC50_STRING as IC50 from NEUT where ANTIBODY = ?) as N
    on N.ID = S.ID and G.ID = S.ID order by S.ID;
    ''', (antibody,))

    # make sure the records are unique
    ab_records = list()
    for row in curr:
        ab_record = ABRecord(*row)
        if ab_record.id in [abr.id for abr in ab_records]:
            ab_record.id += '-1'
        ab_records.append(ab_record)

    conn.close()

    return ab_records


def generate_alignment_from_SeqRecords(filename, seq_records):
    ab_fasta_filename = mkstemp()[1]
    hmm_filename = mkstemp()[1]
    sto_filename = mkstemp()[1]
    finished = False

    try:
        # get the FASTA format file so we can HMMER it
        fafh = open(ab_fasta_filename, 'w')
        hxb2fh = open(OPTIONS.HXB2_FASTA, 'rU')

        # grab the HXB2 Reference Sequence
        hxb2_record = SeqIO.parse(hxb2fh, 'fasta')
        seq_records.extend(hxb2_record)

        # type errors here?
        try:
            SeqIO.write(seq_records, fafh, 'fasta')
        except TypeError, e:
            print >> sys.stderr, seq_records
            raise e

        # close the handles in this order because BioPython wiki does so
        fafh.close()
        hxb2fh.close()

        # make the tempfiles for the alignments, and close them out after we use them
        SeqIO.convert(OPTIONS.HXB2_FASTA, 'fasta', sto_filename, 'stockholm')

        hmmer_build_args = [OPTIONS.HMMER_BUILD_BIN, '--dna' if OPTIONS.DNA else '--amino', \
                           hmm_filename, sto_filename]
        hmmer_align_args = [OPTIONS.HMMER_ALIGN_BIN, '--dna' if OPTIONS.DNA else '--amino', \
                           '--outformat', 'Pfam', '-o', sto_filename, hmm_filename, ab_fasta_filename]

        print >> sys.stderr, 'Aligning %d sequences with HMMER:' % (len(seq_records)+1),

        hmmer = Hmmer(OPTIONS.HMMER_ALIGN_BIN, OPTIONS.HMMER_BUILD_BIN)

        for i in xrange(0, OPTIONS.HMMER_ITER):
            print >> sys.stderr, '%d,' % i,
            hmmer.build(hmm_filename, sto_filename)
            hmmer.align(hmm_filename, ab_fasta_filename, output=sto_filename, alphabet=DNA if OPTIONS.DNA else AMINO, outformat=PFAM)

        # rename the final alignment to its destination
        print >> sys.stderr, 'done, output moved to: %s' % filename
        finished = True

    finally:
        # cleanup these files
        if finished:
            rename(sto_filename, filename)
        else:
            remove(sto_filename)
        remove(ab_fasta_filename)
        remove(hmm_filename)


def clamp(x):
    if x < 0.:
        return 0.
    if x > 1.:
        return 1.
    return x


class SimulatedEpitope(object):

    def __init__(self, positions, position_names, alphabet, kernel_func=None):
        self.positions = positions
        self.alphabet = alphabet
        self.names = position_names
        # default to a uniform linear kernel
        if kernel_func is None:
            self.kernel_func = lambda x, n: clamp(1. * x / len(positions) + n)

    def __str__(self):
        return '\n'.join(sorted(self.names, key = lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x))))

    def evaluate(self, seq, noise=0., proportion=-1):
        total = 0
        # sanitize the sequence, so that it fits in our alphabet
        seq = sanitize_seq(seq, self.alphabet)
        for k, v in self.positions.items():
            if self.alphabet[seq[k]] == self.alphabet[v]:
                total += 1
        # this thing should now produce 50/50 splits no matter what
        # if the positions are more mutated than the base rate, then 25 (resistant)
        # else 1 (susceptible)
        if proportion < 0:
            ret = self.kernel_func(total, noise)
        else:
            ret = 1. if proportion < self.kernel_func(total, noise) else 25.
        # 12.5 * pow(1. - self.kernel_func(total), 0.8) # 0.8 for a 0.1 mutation rate, 2.667 for a 50%/50% split
        # ret = abs(2. * log10(self.kernel_func(total) + 0.00018) / log10(proportion))
        return ret


def random_column_subset(size, columns):
    col_subset = []
    assert(len(columns) > 0)
    while len(col_subset) < size:
        c_ = columns[randint(0, len(columns) - 1)]
        if c_ not in col_subset:
            col_subset.append(c_)
    return col_subset


def generate_random_epitope(seq_table, column_names, size, kernel_func=None):
    if not seq_table.loaded:
        seq_table.fill_columns()

    alphabet, alphabet_names = fetch_alphabet_dict()
    alphabet_len = len(set(alphabet.values()))
    positions = {}

    if size > 0.2 * seq_table.num_columns:
        print >> sys.stderr, 'ERROR: We do not suggest or support simulated epitopes larger than 20% of your alignment length'
        sys.exit(-1)

    # generate a approximately uniformly random epitope from the available nucleotides at each position (ensure that existing sequences match the epitope)
    while len(positions) < size:
        # only grab as many new positions as we need (size - len(positions))
        new_positions = random_column_subset(size - len(positions), seq_table.columns.keys())
        for i in new_positions:
            # if we already have that position, delete it and skip
            if i in positions:
                continue
            # this is some trickery to avoid ambiguous positions and spaces
            vals = [(sanitize_seq(x[0], alphabet), x[1]) for x in seq_table.columns[i].counts().items() if sanitize_seq(x[0], alphabet) not in ('X', '-')]
            # if only Xs, then ignore this position
            if len(vals) == 0:
                continue
            count = sum([x[1] for x in vals])
            vals = sorted([(x[0], 1.0 * x[1] / count) for x in vals], key = itemgetter(1), reverse = True)
            # don't bother with this uniform bullshit, just assume the most common is the epitope
            positions[i] = vals[0][0]
            # find the value whose contribution to the cdf bounds our uniformly random value (r_) , then stop
            # r_ = random()
            # cdf = 0.
            # for j in vals:
                # cdf += j[1]
                # if cdf > r_:
                    # positions[i] = j[0]
                    # break

    # this formula should generate the correct position names for the epitope
    position_names = [column_names[k * alphabet_len + alphabet_names.index(v)] for k, v in positions.items()]

    epi_def = SimulatedEpitope(positions, position_names, alphabet, kernel_func)

    return epi_def


def fetch_alphabet_dict():
    if OPTIONS.STANFEL:
        names = []
        for v in set(_STANFEL.values()):
            names.append('[%s]' % ''.join(sorted([k if k != '-' else '' for (k, v_) in _STANFEL.items() if v == v_])))
        return (_STANFEL, names)
    elif OPTIONS.DNA:
        return (dict([(DNA_ALPHABET[i], i) for i in xrange(0, len(DNA_ALPHABET))]), \
                [DNA_ALPHABET[i] if AMINO_ALPHABET[i] != '-' else '[]' for i in xrange(0, len(DNA_ALPHABET))])
    else:
        return (dict([(AMINO_ALPHABET[i], i) for i in xrange(0, len(AMINO_ALPHABET))]), \
                [AMINO_ALPHABET[i] if AMINO_ALPHABET[i] != '-' else '[]' for i in xrange(0, len(AMINO_ALPHABET))])


def sanitize_seq(seq, alphabet):
    assert(len(SPACE) > 0 and len(seq) > 0 and len(alphabet.keys()) > 0)
    try:
        seq = str(seq)
        seq = seq.upper()
        seq = sub(r'[%s]' % SPACE, '-', seq)
        seq = sub(r'[^%s]' % ''.join(alphabet.keys()), 'X', seq)
    except TypeError, e:
        print >> sys.stderr, 'ERROR: something amiss with things:'
        print >> sys.stderr, 'SPACE =', SPACE
        print >> sys.stderr, 'seq =', seq
        print >> sys.stderr, 'alphabet =', alphabet
        raise e
    return seq


def binarize_row(row):
    alphabet = fetch_alphabet_dict()[0]
    alphabet_len = len(set(alphabet.values()))
    row = sanitize_seq(str(row), alphabet)
    ret = []
    for p in row:
        ret.extend([1 if i == alphabet[p] else 0 for i in xrange(0, alphabet_len)])
    return ret


def generate_feature_names(seq_table):
    alphabet, alphabet_names = fetch_alphabet_dict()

    # make the feature names
    hxb2_seq = None
    for r in seq_table.rows():
        if is_HXB2(r.id):
            hxb2_seq = r.seq

    assert(hxb2_seq is not None)

    # convention is E215 (Glutamic Acid at 215) or 465a (first insertion after 465)
    names = []
    c = 0
    ins = 0
    for p in hxb2_seq:
        if p not in SPACE:
            c += 1
            ins = 0
        else:
            ins += 1
        for v in set(alphabet.values()):
            insert = base_26_to_alph(base_10_to_n(ins, BASE_ALPH))
            names.append('%s%d%s%s' % ('' if insert != '' else p.upper(), c, insert, alphabet_names[v]))

    return names


def compute_relevant_features_and_names(seq_table, names):
    if not seq_table.loaded:
        seq_table.fill_columns()

    columns = range(0, seq_table.num_columns)
    alphabet, alphabet_names = fetch_alphabet_dict()
    alphabet_len = len(set(alphabet.values()))
    delete_cols = list()

    max_counts = list()
    min_counts = list()

    # remove overly-conserved or overly-random columns
    for i in sorted(seq_table.columns.keys()):
        j = alphabet_len * i
        alph_counts = seq_table.columns[i].counts()
        count = sum(alph_counts.values())
        max_count = 1. * max(alph_counts.values()) / count
        min_count = 1. * min(alph_counts.values()) / count
        gap_count = 1. * alph_counts['-'] / count if '-' in alph_counts else 0.
        if max_count > OPTIONS.MAX_CONSERVATION or \
           min_count > OPTIONS.MIN_CONSERVATION or \
           gap_count > OPTIONS.MAX_GAP_RATIO:
            delete_cols.extend(range(j, j+alphabet_len))
            continue
        max_counts.append((i, OPTIONS.MAX_CONSERVATION - max_count))
        min_counts.append((i, OPTIONS.MIN_CONSERVATION - min_count))
        # for each unique assignment
        for v in set(alphabet.values()):
            c = False
            # for each value party to that assignment
            for (k, v_) in alphabet.items():
                if v != v_:
                    continue
                # if we've collected an count for that value
                if k in seq_table.columns[i].counts():
                    c = True
                    break
            # if we've not collected a count, delete it
            if not c:
                delete_cols.append(j+v)

    # TODO: remove overly conserved or overly-random columns in class subgroups

    # I love list comprehensions
    columns = [i for i in xrange(0, seq_table.num_columns * alphabet_len) if i not in delete_cols]

    # trim columns in excess of _MRMR_MAX_VARS or else everything 'splodes
#     if len(columns) * seq_table.num_rows > _MRMR_MAX_VARS:
#         remainder = ceil( (len(columns) * seq_table.num_rows - _MRMR_MAX_VARS) / seq_table.num_rows )
#         print >> sys.stderr, 'WARNING: having to trim %i excess columns, this may take a minute' % remainder
#         max_counts = sorted(max_counts, key=itemgetter(1))
#         max_counts = sorted(min_counts, key=itemgetter(1))
#         while(len(columns) * seq_table.num_rows > _MRMR_MAX_VARS):
#             if max_counts[0][1] < min_counts[0][1]:
#                 j = max_counts.pop(0)[0] * alphabet_len
#                 columns = [i for i in columns if i not in range(j, j+alphabet_len)]
#             else:
#                 j = min_counts.pop(0)[0] * alphabet_len
#                 columns = [i for i in columns if i not in range(j, j+alphabet_len)]
#
#     try:
#         assert(len(columns) <= _MRMR_MAX_VARS)
#     except AssertionError, e:
#         print len(columns), len(delete_cols)
#         raise e

    columns = sorted(columns)

    # trimmed because they're only the ones for columns, not because of the strip
    trimmed_names = [names[i] for i in columns]

    if len(OPTIONS.FILTER) != 0:
        delete_cols = []
        for i in xrange(0, len(columns)):
            m = match(r'^([A-Z]\d+)', trimmed_names[i])
            if m:
                if m.group(1) not in OPTIONS.FILTER:
                    delete_cols.append(i)
        for i in sorted(delete_cols, reverse=True):
            del columns[i]
            del trimmed_names[i]

    return (columns, trimmed_names)


def generate_relevant_data(seq_table, feature_names, feature_idxs, target, proportion=None):
    data = SmlData(feature_names)

    if OPTIONS.SIM in (_RAND_SEQ, _RAND_TARGET) and proportion is None:
        print >> sys.stderr, 'ERROR: no proportion defined for simulation, aborting!'
        sys.exit(-1)

    neg = 0
    pos = 0

    for row in seq_table.rows():
        # everything before the continue takes care of the | and : separators
        if is_HXB2(row.id):
            continue
        if len(OPTIONS.SUBTYPES) != 0:
            subtypes = row.id.split('|')[1].upper()
            if subtypes == '' or len([i for i in subtypes if i in OPTIONS.SUBTYPES]) <= 0:
                print >> sys.stderr, 'ERROR: We\'re supposed to have already masked out unwanted subtypes: %s' % row.id
                raise ValueError
        expanded_row = binarize_row(row.seq)
        feats = dict([f for f in zip(range(0, len(feature_idxs)), [expanded_row[i] for i in feature_idxs]) if f[1] != 0])
        if OPTIONS.SIM in (_RAND_SEQ, _RAND_TARGET):
            class_ = random()
            # if the proportion is 1, then 100% of class_ should be 1, if proportion is .3, then 30%, and so on
            class_ = 1 if class_ < proportion else 0
        else:
            class_ = id_to_class(row.id, target)
        data.add(class_, feats)
        if class_:
            pos += 1
        else:
            neg += 1

    # silence this during testing
    if not OPTIONS.TEST and 0:
        print >> sys.stdout, 'Ratio + to - : %d to %d' % (pos, neg)

    return data


def collect_seq_table_and_feature_names(ab_alignment_filename, seq_records):

    if not exists(ab_alignment_filename) and OPTIONS.SIM != _RAND_DUMB:
        generate_alignment_from_SeqRecords(ab_alignment_filename, seq_records)
    elif OPTIONS.SIM == _RAND_DUMB:
        # stupid workaround because Bio.SeqIO objects suck at regular Stockholm output
        ab_alignment_fh = open(ab_alignment_filename, 'w')
        SeqIO.write(seq_records, ab_alignment_fh, 'stockholm')
        ab_alignment_fh.close()

    skip_func = None
    if len(OPTIONS.SUBTYPES) != 0:
        skip_func = lambda x: len([i for i in x.split('|')[1] if i in tuple(OPTIONS.SUBTYPES)]) <= 0

    seq_table = SeqTable(ab_alignment_filename, DNA_ALPHABET if OPTIONS.DNA else AMINO_ALPHABET, is_HXB2, skip_func)

    all_feature_names = generate_feature_names(seq_table)

    # set the right number of CV_FOLDS for Leave-One-Out-Crossvalidation
    if OPTIONS.LOOCV:
        OPTIONS.CV_FOLDS = len([r for r in seq_table.rows() if not is_HXB2(r.id)])

    # parition the table into CV_FOLDS groups for cross-validation
    seq_table.partition(OPTIONS.CV_FOLDS)

    return seq_table, all_feature_names


def assign_class_by_percentile(seq_table, epi_def):

    vals = []
    for row in seq_table.rows():
        if is_HXB2(row.id):
            continue
        vals.append(epi_def.evaluate(row.seq, get_noise(row.id)))
    proportion = percentile(vals, OPTIONS.SIM_EPI_PERCENTILE)

    for row in seq_table.rows():
        if is_HXB2(row.id):
            continue
        row.id = '|||%.3f' % epi_def.evaluate(row.seq, get_noise(row.id), proportion)

    return


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
    valid_antibodies = sorted(get_valid_antibodies_from_db(), key = lambda x: x.strip())
    if antibody not in valid_antibodies:
        if ' ' + antibody not in valid_antibodies:
            option_parser.error('%s not in the list of permitted antibodies: \n  %s' % (antibody, '\n  '.join([ab.strip() for ab in valid_antibodies])))
        else:
            antibody = ' ' + antibody

    # validate the subtype option
    valid_subtypes = sorted(get_valid_subtypes_from_db(), key = lambda x: x.strip().upper())
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

    ab_basename = '%s%s_%s' % (antibody, '_randseq' if OPTIONS.SIM == _RAND_SEQ else '', 'dna' if OPTIONS.DNA else 'amino')

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    ab_records = collect_ABRecords_from_db(antibody)

    # do not generate the sequences here for the random sequence or random epitope simulations
    ab_alignment_filename = None
    if OPTIONS.SIM not in (_RAND_DUMB, _RAND_EPI, _RAND_SEQ):
        ab_alignment_filename = '%s_%s_%s.sto' % (ab_basename, splitext(basename(OPTIONS.NEUT_SQLITE3_DB))[0], __version__)

        # generate an alignment using HMMER if it doesn't already exist
        seq_records = [r.to_SeqRecord() for r in ab_records]
        seq_table, all_feature_names = collect_seq_table_and_feature_names(ab_alignment_filename, seq_records)
    else:
        if OPTIONS.SIM_EPI_N is None:
            OPTIONS.SIM_EPI_N = len(ab_records)

    # fetch the alphabet, we'll probably need it later
    alphabet = fetch_alphabet_dict()[0]

    # compute features
    for target in OPTIONS.TARGETS:
        results = None

        if OPTIONS.SIM not in (_RAND_DUMB, _RAND_EPI, _RAND_SEQ):
            seq_table.unmask()

        # simulations, ho!
        for i in xrange(0, OPTIONS.SIM_RUNS if OPTIONS.SIM != '' else 1):

            # here is where the sequences must be generated for the random sequence and random epitope simulations
            if OPTIONS.SIM in (_RAND_DUMB, _RAND_EPI, _RAND_SEQ):
                ab_alignment_filename = '%s_%s_%d_%s.sto' % (ab_basename, OPTIONS.SIM, i + 1, __version__)
                if OPTIONS.SIM == _RAND_DUMB:
                    hxb2fh = open(OPTIONS.HXB2_FASTA)
                    for record in SeqIO.parse(hxb2fh, 'fasta'):
                        hxb2_record = record
                    seq_records = DumbRandomSequences(str(hxb2_record.seq), idfmt='%s|||%.3g', N=OPTIONS.SIM_EPI_N, noise=OPTIONS.SIM_EPI_NOISE, rate=OPTIONS.SIM_EPI_MUT_RATE, alphabet=alphabet)
                    seq_records.append(hxb2_record)
                    hxb2fh.close()
                else:
                    seq_records = MarkovRandomSequences(_RAND_SEQ_STOCKHOLM, idfmt='%s|||', N=len(ab_records), noise=OPTIONS.SIM_EPI_NOISE, rate=OPTIONS.SIM_EPI_MUT_RATE, alphabet=alphabet)
                seq_table, all_feature_names = collect_seq_table_and_feature_names(ab_alignment_filename, seq_records)

            if OPTIONS.SIM in (_RAND_DUMB, _RAND_EPI):
                epi_def = generate_random_epitope(seq_table, all_feature_names, OPTIONS.SIM_EPI_SIZE)
                print >> sys.stdout, '********************* SIMULATED EPITOPE DESCRIPTION (%d) *********************\n' % OPTIONS.SIM_EPI_SIZE
                print >> sys.stdout, '%s\n' % str(epi_def)
                # assign the appropriate class by percentile
                assign_class_by_percentile(seq_table, epi_def)
            else:
                epi_def = None

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
                selectorcls=Mrmr,
                folds=OPTIONS.CV_FOLDS,
                cv={},
                optstat=optstat,
                gs={ 'C': C_range },
                fs={ 'num_features': OPTIONS.NUM_FEATURES, 'method': Mrmr.MAXREL if OPTIONS.MAXREL else Mrmr.MID if OPTIONS.MRMR_METHOD == 'MID' else Mrmr.MIQ }
            )

            # make sure the whole thing is unmasked for the nestedcrossvalidator
            seq_table.unmask()
            feature_idxs, feature_names = compute_relevant_features_and_names(seq_table, all_feature_names)
            smldata = generate_relevant_data(seq_table, feature_names, feature_idxs, target)
            x, y = smldata.tondarrays()

            results = crossvalidator.crossvalidate(x, y, cv={}, extra=lambda x: { 'features': x.features(), 'weights': x.classifier.weights() })

            # remove the alignment
            if OPTIONS.SIM in (_RAND_DUMB, _RAND_EPI):
                remove(ab_alignment_filename)

        # print stats and results:
        print >> sys.stdout, '********************* REPORT FOR ANTIBODY %s IC50 %s *********************' % \
          (antibody, '< %d' % OPTIONS.IC50LT if target == 'lt' else '> %d' % OPTIONS.IC50GT)

        if OPTIONS.SIM not in (_RAND_SEQ, _RAND_TARGET, _RAND_EPI):
            fmt = ('', OPTIONS.CV_FOLDS, '')
        else:
            fmt = ('%d-run ' % OPTIONS.SIM_RUNS, OPTIONS.CV_FOLDS, ' per run')

        print >> sys.stdout, '\n%sLSVM performance statistics (%d-fold CV%s):' % fmt

        statsdict = results['stats'].todict()

        # remove minstat 'cause we don't want it here.. 
        if 'Minstat' in statsdict:
            del statsdict['Minstat']

        # convert to percentage from [0, 1]
        for k, v in statsdict.items():
            statsdict[k] = v * 100.

        stat_len = max([len(k) for k in statsdict.keys()])
        mean_len = max([len('%.2f' % v.mu) for v in statsdict.values()])
        std_len = max([len('%.2f' % sqrt(v.sigma)) for v in statsdict.values()]) 
        fmt = u'%%%d.2f \xb1 %%%d.2f%%%%' % (mean_len, std_len) 
        for k, v in sorted(statsdict.items(), key=itemgetter(0)):
            print >> sys.stdout, u'  %-*s = %s' % (stat_len, k, v.sprintf(fmt)) 

        features = set()
        featureweights = []
        for i in xrange(len(results['extra'])):
            assert(len(results['extra'][i]['features']) >= len(results['extra'][i]['weights']))
            for j in xrange(len(results['extra'][i]['features'])):
                w = results['extra'][i]['weights'][j] if j < len(results['extra'][i]['weights']) else 0.
                featureweights.append((results['extra'][i]['features'][j], w)) 

        weightsdict = {}
        for featureidx in set([fw[0] for fw in featureweights]):
            weights = [int(copysign(1, fw[1])) for fw in featureweights if fw[0] == featureidx]
            val = NormalValue(int, weights)
            if (abs(val.mu) < 0.0001 and val.sigma == 0.) or len(val) == 1:
                continue
            weightsdict[feature_names[featureidx]] = val

        print >> sys.stdout, '\nSignificant positions (top %d):' % (len(weightsdict))

        if len(weightsdict) > 0:
            name_len = max([len(k) for k in weightsdict.keys()])
            mean_len = max([len('% .1f' % v.mu) for v in weightsdict.values()])
            std_len = max([len('%.1f' % sqrt(v.sigma)) for v in weightsdict.values()])
            fmt = u'%% %d.1f \xb1 %%%d.1f' % (mean_len, std_len)
            N_len = max([len('%d' % len(v)) for v in weightsdict.values()])
            for k, v in sorted(weightsdict.items(), key = lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x[0]))):
                print >> sys.stdout, u'  %-*s  %s  (N = %*d)' % (name_len, k, v.sprintf(fmt), N_len, len(v))

        print >> sys.stdout, '\n'

        # TODO: perform patch analysis

    return 0


if __name__ == '__main__':
    sys.exit(main())
