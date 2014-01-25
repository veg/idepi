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

from __future__ import division, print_function

import sys

from optparse import OptionParser
from os import close, remove
from os.path import abspath, join, split
from random import seed
from re import search
from tempfile import mkstemp

import numpy as np

from Bio import AlignIO

from idepi.alphabet import Alphabet
from idepi.databuilder import DataBuilder
from idepi.filter import naivefilter

from idepi import (
    Labeler,
    Regressor,
    alignment_identify_ref,
    extract_feature_weights,
    cv_results_to_output,
    generate_alignment,
    input_data,
    is_HXB2,
    regressor_classes,
    pretty_fmt_results,
    seqrecord_get_values,
    set_util_params,
    __file__ as _idepi_file,
    __version__ as _idepi_version
)

from pyxval import ContinuousPerfStats, CrossValidator

__VERSION__ = _idepi_version

_IDEPI_PATH       = split(abspath(_idepi_file))[0]
_HXB2_DNA_FASTA   = join(_IDEPI_PATH, 'data', 'hxb2_dna.fa')
_HXB2_AMINO_FASTA = join(_IDEPI_PATH, 'data', 'hxb2_pep.fa')

_PHYLOFILTER_BATCHFILE = join(_IDEPI_PATH, 'data', 'hyphy', 'CorrectForPhylogeny.bf')

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

_TEST_Y = np.array([1.,21.,50.,0.5])

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


def optparse_data(option, _, value, parser):
    setattr(parser.values, option.dest, input_data(value))


def setup_option_parser():

    parser = OptionParser(usage = '%prog [options] ANTIBODY')

    #                 option            action='store'        callback                 type             nargs=1  dest
    parser.add_option('--hmmalign',                                                    type='string',            dest='HMMER_ALIGN_BIN')
    parser.add_option('--hmmbuild',                                                    type='string',            dest='HMMER_BUILD_BIN')
    parser.add_option('--hmmiter',                                                     type='int',               dest='HMMER_ITER')
    parser.add_option('--method',                                                      type='string',            dest='REGRESSOR_METHOD')
    parser.add_option('--filter',       action='callback',    callback=optparse_csv,   type='string',            dest='FILTER')
    parser.add_option('--clonal',       action='store_true',                                                     dest='CLONAL')
    parser.add_option('--numfeats',                                                    type='int',               dest='NUM_FEATURES')
    parser.add_option('--subtypes',     action='callback',    callback=optparse_csv,   type='string',            dest='SUBTYPES')
    parser.add_option('--weighting',    action='store_true',                                                     dest='WEIGHTING')
    parser.add_option('--amino',        action='store_true',                                                     dest='AMINO')
    parser.add_option('--dna',          action='store_true',                                                     dest='DNA')
    parser.add_option('--stanfel',      action='store_true',                                                     dest='STANFEL')
    parser.add_option('--cv',                                                          type='int',               dest='CV_FOLDS')
    parser.add_option('--loocv',        action='store_true',                                                     dest='LOOCV')
    parser.add_option('--maxcon',                                                      type='float',             dest='MAX_CONSERVATION')
    parser.add_option('--maxgap',                                                      type='float',             dest='MAX_GAP_RATIO')
    parser.add_option('--mincon',                                                      type='float',             dest='MIN_CONSERVATION')
    parser.add_option('--neuts',        action='callback',    callback=optparse_data,  type='string',            dest='DATA')
    parser.add_option('--hxb2',                                                        type='string',            dest='REFSEQ_FASTA')
    parser.add_option('--ids',          action='callback',    callback=optparse_csv,   type='string',            dest='HXB2_IDS')
    parser.add_option('--test',         action='store_true',                                                     dest='TEST')
    parser.add_option('--seed',                                                        type='int',               dest='RAND_SEED')
    parser.add_option('--phylofilt',    action='store_true',                                                     dest='PHYLOFILTER')
    parser.add_option('--logspace',     action='store_true',                                                     dest='LOGSPACE')

    parser.set_defaults(HMMER_ALIGN_BIN    = 'hmmalign')
    parser.set_defaults(HMMER_BUILD_BIN    = 'hmmbuild')
    parser.set_defaults(HMMER_ITER         = 8)
    parser.set_defaults(REGRESSOR_METHOD   = 'ridgelar')
    parser.set_defaults(FILTER             = [])
    parser.set_defaults(CLONAL             = False)
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
    parser.set_defaults(DATA               = input_data(join(_IDEPI_PATH, 'data', 'allneuts.sqlite3')))
    parser.set_defaults(REFSEQ_FASTA       = _HXB2_AMINO_FASTA)
    parser.set_defaults(HXB2_IDS           = ['9629357', '9629363'])
    parser.set_defaults(RAND_SEED          = 42) # make the behavior deterministic for now
    parser.set_defaults(PHYLOFILTER        = False)
    parser.set_defaults(LOGSPACE           = False)

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

            # test mRMR and LSVM file generation
            ylabeler = Labeler(
                seqrecord_get_values,
                lambda seq: is_HXB2(seq) or False, # TODO: again filtration function
            )
            alignment, y, ic50gt = ylabeler(alignment)

            filter = naivefilter(
                OPTIONS.MAX_CONSERVATION,
                OPTIONS.MIN_CONSERVATION,
                OPTIONS.MAX_GAP_RATIO
            )
            refidx = alignment_identify_ref(alignment, is_HXB2)
            builder = DataBuilder(
                alignment,
                alph,
                refidx,
                filter
            )
            x = builder(alignment, refidx)
            colnames = builder.labels

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

            assert(np.all(_TEST_Y == y))

            # TODO: generate and test the regressor data generation
            # print y, "\n", x

    finally:
        remove(sto_filename)

    print('ALL TESTS PASS', file=sys.stderr)


def fix_hxb2_fasta():
    '''If DNA mode was selected but the AMINO reference sequence is still in place, fix it'''
    if OPTIONS.DNA == True and OPTIONS.REFSEQ_FASTA == _HXB2_AMINO_FASTA:
        OPTIONS.REFSEQ_FASTA = _HXB2_DNA_FASTA


def main(argv=sys.argv):
    global OPTIONS

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
    cvopts = {}
    if OPTIONS.REGRESSOR_METHOD in regressor_classes:
        cvopts['regressorcls'] = regressor_classes[OPTIONS.REGRESSOR_METHOD]
    else:
        option_parser.error('%s not in the list of available regression methods: \n  %s' % (OPTIONS.REGRESSOR_METHOD,
            '\n  '.join(regressor_classes.keys())))

    if search(r'(?:lar|lasso)$', OPTIONS.REGRESSOR_METHOD):
        if OPTIONS.NUM_FEATURES < 0:
            OPTIONS.NUM_FEATURES = _DEFAULT_NUM_FEATURES
        cvopts['m'] = OPTIONS.NUM_FEATURES
    elif OPTIONS.NUM_FEATURES > 0:
        option_parser.error('--numfeats is a useless parameter for regression method `%s\'' % OPTIONS.REGRESSOR_METHOD)

    cvopts['logspace'] = OPTIONS.LOGSPACE

    # validate the antibody argument, currently a hack exists to make PG9/PG16 work
    # TODO: Fix pg9/16 hax
    antibody = args[1].strip()
    valid_antibodies = sorted(OPTIONS.DATA.antibodies, key=lambda x: x.strip())
    if antibody not in valid_antibodies:
        if ' ' + antibody not in valid_antibodies:
            option_parser.error('%s not in the list of permitted antibodies: \n  %s' % (antibody, '\n  '.join([ab.strip() for ab in valid_antibodies])))
        else:
            antibody = ' ' + antibody

    # validate the subtype option
    valid_subtypes = sorted(OPTIONS.DATA.subtypes, key=lambda x: x.strip().upper())
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

    # fetch the alphabet, we'll probably need it later
    alph = Alphabet(mode=Alphabet.STANFEL if OPTIONS.STANFEL else Alphabet.DNA if OPTIONS.DNA else Alphabet.AMINO)

    ab_basename = ''.join((
        antibody,
        '_dna' if OPTIONS.DNA else '_amino',
        '_clonal' if OPTIONS.CLONAL else ''
    ))
    alignment_basename = '_'.join((
        ab_basename,
        OPTIONS.DATA.basename_root,
        __VERSION__
    ))

    # grab the relevant antibody from the SQLITE3 data
    # format as SeqRecord so we can output as FASTA
    # and generate an alignment using HMMER if it doesn't already exist
    seqrecords, clonal = OPTIONS.DATA.seqrecords(antibody, OPTIONS.CLONAL, OPTIONS.DNA)

    # if clonal isn't supported, fallback to default
    if clonal != OPTIONS.CLONAL:
        ab_basename = ''.join(ab_basename.rsplit('_clonal', 1))
        alignment_basename = ''.join(alignment_basename.rsplit('_clonal', 1))

    sto_filename = alignment_basename + '.sto'

    alignment = generate_alignment(seqrecords, sto_filename, is_refidx, OPTIONS)[0]

    ylabeler = Labeler(
        seqrecord_get_values,
        lambda seq: is_HXB2(seq) or False, # TODO: again filtration function
    )
    alignment, y, ic50gt = ylabeler(alignment)

    filter = naivefilter(
        OPTIONS.MAX_CONSERVATION,
        OPTIONS.MIN_CONSERVATION,
        OPTIONS.MAX_GAP_RATIO,
    )
    refidx = alignment_identify_ref(alignment, is_HXB2)
    builder = DataBuilder(
        alignment,
        alph,
        refidx,
        filter
    )
    x = builder(alignment, refidx)
    colnames = builder.labels

    crossvalidator = CrossValidator(
        classifier_cls=Regressor,
        folds=OPTIONS.CV_FOLDS,
        classifier_kwargs=cvopts,
        scorer_cls=ContinuousPerfStats,
        scorer_kwargs={}
    )

    results = crossvalidator.crossvalidate(x, y, classifier_kwargs={}, extra=extract_feature_weights)

    ret = cv_results_to_output(results, colnames)

    print(pretty_fmt_results(ret))

#     mean_len = max([len('%.3f' % v.mu) for v in avg_stats.values()])
#     std_len = max([len('%.3f' % v.sigma) for v in avg_stats.values()])
#     std_len = int(log10(max([1.] + [v.sigma for v in avg_stats.values()]))) + 5
#     for k, v in sorted(avg_stats.items(), key = lambda x: x[0][0]):
#         v_str = u'= %*.3f \xb1 %*.3f' % (mean_len, v.mu, std_len, v.sigma)
#         print(u'  %s%s' % (k, v_str))
#
#     for k, v in avg_weights.items():
#         if abs(v.mu) < 0.0001 and v.sigma == 0.:
#             del avg_weights[k]
#
#     print('\nSignificant positions (top %d):' % (len(avg_weights)))
#
#     if len(avg_weights) > 0:
#         name_len = max(len(k) for k in avg_weights.keys())
#         mean_len = max(len('% .1f' % v.mu) for v in avg_weights.values())
#         std_len = max(len('%.1f' % v.sigma) for v in avg_weights.values())
#         N_len = max(len('%d' % len(v.values)) for v in avg_weights.values())
#         for k, v in sorted(avg_weights.items(), key=lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x[0]))):
#             print(u'  %-*s  % *.1f \xb1 %*.1f (N = %*d)' % (name_len, k, mean_len, v.mu, std_len, v.sigma, N_len, len(v.values)))
#
#     print('\n')

    return 0


if __name__ == '__main__':
    sys.exit(main())
