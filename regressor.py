#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# regressor.py :: computes a regression model for an nAb in IDEPI's database, 
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
from math import ceil, log10
from operator import itemgetter
from optparse import OptionParser
from os import remove, rename
from os.path import basename, dirname, exists, join, realpath, splitext
from random import gauss, randint, random, seed
from re import match, search, sub
from subprocess import Popen, PIPE
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

from idepi import REGRESSOR_METHODS, Regressor
from idepi import OrfList
from idepi import DumbRandomSequences, MarkovRandomSequences
from idepi import AMINO_ALPHABET, DNA_ALPHABET, SeqTable, SPACE
from idepi import SmlData
from idepi import LinearSvm, NormalValue
from idepi import get_noise, set_util_params, is_HXB2, id_to_real, percentile, base_10_to_n, base_26_to_alph, BASE_ALPH

from numpy import mean, median, std

__version__ = 0.4

_WORKING_DIR      = dirname(realpath(__file__)) # '/Users/Lance/Pond'
_HXB2_DNA_FASTA   = join(_WORKING_DIR, 'res', 'hxb2_dna.fa')
_HXB2_AMINO_FASTA = join(_WORKING_DIR, 'res', 'hxb2_pep.fa')

_DEFAULT_NUM_FEATURES = 10

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

    #                 option           action = 'store'        callback                    type             nargs = 1  dest
    parser.add_option('--hmmalign',                                                         type = 'string',            dest = 'HMMER_ALIGN_BIN')
    parser.add_option('--hmmbuild',                                                         type = 'string',            dest = 'HMMER_BUILD_BIN')
    parser.add_option('--hmmiter',                                                          type = 'int',               dest = 'HMMER_ITER')
    parser.add_option('--method',                                                           type = 'string',            dest = 'REGRESSOR_METHOD')
    parser.add_option('--filter',       action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'FILTER')
    parser.add_option('--numfeats',                                                         type = 'int',               dest = 'NUM_FEATURES')
    parser.add_option('--subtypes',     action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'SUBTYPES')
    parser.add_option('--weighting',    action = 'store_true',                                                          dest = 'WEIGHTING')
    parser.add_option('--amino',        action = 'store_true',                                                          dest = 'AMINO')
    parser.add_option('--dna',          action = 'store_true',                                                          dest = 'DNA')
    parser.add_option('--stanfel',      action = 'store_true',                                                          dest = 'STANFEL')
    parser.add_option('--cv',                                                               type = 'int',               dest = 'CV_FOLDS')
    parser.add_option('--loocv',        action = 'store_true',                                                          dest = 'LOOCV')
    parser.add_option('--maxcon',                                                           type = 'float',             dest = 'MAX_CONSERVATION')
    parser.add_option('--maxgap',                                                           type = 'float',             dest = 'MAX_GAP_RATIO')
    parser.add_option('--mincon',                                                           type = 'float',             dest = 'MIN_CONSERVATION')
    parser.add_option('--neuts',                                                            type = 'string',            dest = 'NEUT_SQLITE3_DB')
    parser.add_option('--hxb2',                                                             type = 'string',            dest = 'HXB2_FASTA')
    parser.add_option('--ids',          action = 'callback',    callback = optparse_csv,    type = 'string',            dest = 'HXB2_IDS')
    parser.add_option('--test',         action = 'store_true',                                                          dest = 'TEST')
    parser.add_option('--seed',                                                             type = 'int',               dest = 'RAND_SEED')

    parser.set_defaults(HMMER_ALIGN_BIN    = join(_WORKING_DIR, 'contrib', 'hmmer-3.0', 'src', 'hmmalign'))
    parser.set_defaults(HMMER_BUILD_BIN    = join(_WORKING_DIR, 'contrib', 'hmmer-3.0', 'src', 'hmmbuild'))
    parser.set_defaults(HMMER_ITER         = 8)
    parser.set_defaults(REGRESSOR_METHOD   = 'lar')
    parser.set_defaults(FILTER             = [])
    parser.set_defaults(NUM_FEATURES       = -1)
    parser.set_defaults(SUBTYPES           = [])
    parser.set_defaults(WEIGHTING          = False)
    parser.set_defaults(AMINO              = False)
    parser.set_defaults(DNA                = False)
    parser.set_defaults(STANFEL            = False)
    parser.set_defaults(CV_FOLDS           = 5)
    parser.set_defaults(LOOCV              = False)
    parser.set_defaults(MAX_CONSERVATION   = 1. ) # 93.)
    parser.set_defaults(MAX_GAP_RATIO      = 0.1 ) # 93.)
    parser.set_defaults(MIN_CONSERVATION   = 1. ) # 33.)
    parser.set_defaults(NEUT_SQLITE3_DB    = join(_WORKING_DIR, 'res', 'allneuts.sqlite3'))
    parser.set_defaults(HXB2_FASTA         = _HXB2_AMINO_FASTA)
    parser.set_defaults(HXB2_IDS           = ['9629357', '9629363'])
    parser.set_defaults(RAND_SEED          = 42) # make the behavior deterministic for now

    return parser


def run_tests():
    # set these to this so we don't exclude anything (just testing file generation and parsing)
    OPTIONS.NUM_FEATURES = 15 # should be enough, the number is known to be 13
    OPTIONS.MAXREL = False
    OPTIONS.DNA = False
    OPTIONS.MAX_CONSERVATION = 1.0
    OPTIONS.MAX_GAP_RATIO    = 1.0
    OPTIONS.MIN_CONSERVATION = 1.0

    # if we don't do this, DOOMBUNNIES
    set_util_params(OPTIONS.HXB2_IDS)

    fd, sto_filename = mkstemp(); close(fd)

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
                assert(len(feature_names) == len(_TEST_NAMES))
            except AssertionError, e:
                print 'gen:   %s\ntruth: %s' % (feature_names, _TEST_NAMES)
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
            data = generate_relevant_data(seq_table, feature_names, feature_idxs)
    
            # TODO: generate and test the regressor data generation    

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
        self.dna_seq, self.amino_seq = ORFList(seq)[0]
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
    fd, ab_fasta_filename = mkstemp(); close(fd)
    fd, sto_filename = mkstemp(); close(fd)
    fd, hmm_filename = mkstemp(); close(fd)
    finished = False

    try:
        # get the FASTA format file so we can HMMER it
        fafh = open(ab_fasta_filename, 'w')
        hxb2fh = open(OPTIONS.HXB2_FASTA, 'rU')
    
        # grab the HXB2 Reference Sequence
        hxb2_record = SeqIO.parse(hxb2fh, 'fasta')
        seq_records.extend(hxb2_record)
        
        try:
            SeqIO.write(seq_records, fafh, 'fasta')
        except TypeError, e:
            print seq_records
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
    
        for i in xrange(0, OPTIONS.HMMER_ITER):
            print >> sys.stderr, '%d,' % i,
            build_process = Popen(hmmer_build_args, close_fds = True, stderr = PIPE, stdout = PIPE)
            build_process.communicate()
            align_process = Popen(hmmer_align_args, close_fds = True, stderr = PIPE, stdout = PIPE)
            align_process.communicate()
    
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
    # TODO: perform `blurring' based on HMM posterior probabilities for windows of length 3, 5, or 7
    # also investigate greater blurring between amino-acid identities based on some measure of
    # relatedness (BLOSUM?) 
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


def generate_relevant_data(seq_table, feature_names, feature_idxs, target=None, proportion=None):
    data = SMLData(feature_names)

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
        class_ = id_to_real(row.id)
        data.add(class_, feats)

    return data


def collect_seq_table_and_feature_names(ab_alignment_filename, seq_records):

    if not exists(ab_alignment_filename):
        generate_alignment_from_SeqRecords(ab_alignment_filename, seq_records)

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

    if OPTIONS.RAND_SEED is not None:
        seed(OPTIONS.RAND_SEED)

    if len(args) != 2:
        option_parser.error('ANTIBODY is a required argument')

    # check to make sure our mode is exclusive, and set the default (AMINO) if none is set
    if sum([1 for v in (OPTIONS.AMINO, OPTIONS.DNA, OPTIONS.STANFEL) if v]) > 1:
        option_parser.error('options --amino, --dna, and --stanfel are mutually exclusive')
    elif sum([1 for v in (OPTIONS.AMINO, OPTIONS.DNA, OPTIONS.STANFEL) if v]) == 0:
        OPTIONS.AMINO = True

    # validate the regression method
    if OPTIONS.REGRESSOR_METHOD not in regressor_methods.keys():
        option_parser.error('%s not in the list of available regression methods: \n  %s' % (OPTIONS.REGRESSOR_METHOD, '\n  '.join(regressor_methods.keys())))

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

    if len(OPTIONS.FILTER) != 0:
        if OPTIONS.NUM_FEATURES != -1:
            option_parser.error('--filter and --numfeats are incompatible options')
        else:
            OPTIONS.NUM_FEATURES = len(OPTIONS.FILTER)
    else: # len(OPTIONS.FILTER) == 0
        if OPTIONS.NUM_FEATURES == -1:
            OPTIONS.NUM_FEATURES = _DEFAULT_NUM_FEATURES

    # destroy the parser because optparse docs recommend it
    option_parser.destroy()

    # use the default DNA HXB2 Reference seq if we define --dna but don't give a new default HXB2 Reference seq
    fix_hxb2_fasta()

    # set the util params
    set_util_params(OPTIONS.HXB2_IDS)

    ab_basename = '%s_%s' % (antibody, 'dna' if OPTIONS.DNA else 'amino')

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    ab_records = collect_ABRecords_from_db(antibody)

    ab_alignment_filename = '%s_%s_%s.sto' % (ab_basename, splitext(basename(OPTIONS.NEUT_SQLITE3_DB))[0], __version__)
    
    # generate an alignment using HMMER if it doesn't already exist
    seq_records = [r.to_SeqRecord() for r in ab_records]
    seq_table, all_feature_names = collect_seq_table_and_feature_names(ab_alignment_filename, seq_records)

    # fetch the alphabet, we'll probably need it later
    alphabet = fetch_alphabet_dict()[0]

    # compute features
    avg_stats = {}
    avg_weights = {}
 
    seq_table.unmask()
   
    for j in xrange(0, OPTIONS.CV_FOLDS):

        # generate the training data by masking out 1/CV_FOLDS worth
        seq_table.mask(j)
        
        feature_idxs, feature_names = compute_relevant_features_and_names(seq_table, all_feature_names)
        train_data = generate_relevant_data(seq_table, feature_names, feature_idxs)

        # generate the test data by masking out everything but that 1/CV_FOLDS worth
        seq_table.mask([i for i in xrange(0, OPTIONS.CV_FOLDS) if i != j])
        test_data = generate_relevant_data(seq_table, feature_names, feature_idxs)

        regressor_options = {}
        if search(r'(?:lar|lasso)$', OPTIONS.REGRESSOR_METHOD):
            regressor_options['m'] = OPTIONS.NUM_FEATURES
        regressor = Regressor(train_data, method=OPTIONS.REGRESSOR_METHOD, **regressor_options)
        stats = regressor.test(test_data)
        weights = regressor.weights
        
        for k, v in stats.items():
            if k not in avg_stats:
                avg_stats[k] = NormalValue([], k)
            avg_stats[k].append(v)
        
        assert(len(train_data.feature_names) == len(weights))
        for k, v in weights.items():
            if k not in avg_weights:
                avg_weights[k] = NormalValue([], k)
            avg_weights[k].append(v)

    # print stats and results:
    print >> sys.stdout, '********************* REPORT FOR ANTIBODY %s *********************' % \
      (antibody,)

    fmt_vals = (regressor_methods[OPTIONS.REGRESSOR_METHOD], OPTIONS.CV_FOLDS)

    print >> sys.stdout, '\n%s performance statistics (%d-fold CV):' % fmt_vals 
    
    mean_len = max([len('%.3f' % v.mu) for v in avg_stats.values()])
    std_len = max([len('%.3f' % v.sigma) for v in avg_stats.values()])
    std_len = int(log10(max([1.] + [v.sigma for v in avg_stats.values()]))) + 5
    for k, v in sorted(avg_stats.items(), key = lambda x: x[0][0]):
        v_str = u'= %*.3f \xb1 %*.3f' % (mean_len, v.mu, std_len, v.sigma)
        print >> sys.stdout, u'  %s%s' % (k, v_str)

    for k, v in avg_weights.items():
        if abs(v.mu) < 0.0001 and v.sigma == 0.:
            del avg_weights[k]

    print >> sys.stdout, '\nSignificant positions (top %d):' % (len(avg_weights))

    if len(avg_weights) > 0:
        name_len = max([len(k) for k in avg_weights.keys()])
        mean_len = max([len('% .1f' % v.mu) for v in avg_weights.values()])
        std_len = max([len('%.1f' % v.sigma) for v in avg_weights.values()])
        N_len = max([len('%d' % len(v.values)) for v in avg_weights.values()])
        for k, v in sorted(avg_weights.items(), key = lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x[0]))):
            print >> sys.stdout, u'  %-*s  % *.1f \xb1 %*.1f (N = %*d)' % (name_len, k, mean_len, v.mu, std_len, v.sigma, N_len, len(v.values))

    print >> sys.stdout, '\n'

    # TODO: perform patch analysis

    return 0


if __name__ == '__main__':
    sys.exit(main())
