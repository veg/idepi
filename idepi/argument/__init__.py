# _cmdopts.py :: common command options for IDEPI programs
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

# this file is heavily formatted, so skip it for the most part
# flake8: noqa


from __future__ import division, print_function

import logging
import sys

from argparse import ArgumentParser, ArgumentTypeError, FileType
from os.path import isfile, join
from random import seed

from numpy.random import seed as np_seed

from Bio import SeqIO

from BioExt.misc import translate
from BioExt.references import hxb2

from idepi.constants import AminoAlphabet, DNAAlphabet
from idepi.encoder import (
    AminoEncoder,
    DNAEncoder,
    StanfelEncoder
    )
from idepi.datasource import DataSource
from idepi.scorer import Scorer
from idepi.simulation import Simulation
from idepi.util import seqfile_format
from idepi.verifier import VerifyError, Verifier


def RangesType(string):
    msg = "invalid comma-delimited list of numbers '{0:s}'".format(string)
    try:
        r = set()
        for s in string.split(','):
            args = [int(t) for t in s.split(':')]
            if len(args) > 1:
                args[1] += 1  # ranges should be inclusive, for ease of use
                r.update(range(*args))
            else:
                r.update(args)
        return sorted(r)
    except TypeError:
        raise ArgumentTypeError(msg)
    except ValueError:
        raise ArgumentTypeError(msg)


def csvtype(string):
    return string.split(',')


def logtype(string):
    loggers = set(v.lower() for v in string.split(','))
    all = set(['idepi'])
    if 'all' in loggers:
        loggers |= all
    diff = all - loggers
    if len(diff):
        raise ValueError("unknown loggers requested: %s" % ', '.join("'%s'" for l in diff))
    # set the loggers
    if 'idepi' in loggers:
        from idepi.logging import IDEPI_LOGGER
        logging.getLogger(IDEPI_LOGGER).setLevel(logging.DEBUG)
    return loggers


def hmmer_args(parser):
    #                   option        type      dest
    parser.add_argument('--hmmalign', type=str, dest='HMMER_ALIGN_BIN')
    parser.add_argument('--hmmbuild', type=str, dest='HMMER_BUILD_BIN')
    parser.set_defaults(
        HMMER_ALIGN_BIN='hmmalign',
        HMMER_BUILD_BIN='hmmbuild',
        )
    return parser


def featsel_args(parser):
    #                   option             action                type             dest
    parser.add_argument('--numfeats',                            type=RangesType, dest='FEATURE_GRID')
    parser.set_defaults(
        FEATURE_GRID=[10]
        )
    return parser


def mrmr_args(parser):
    group = parser.add_mutually_exclusive_group()
    #                   option             action                const              type        dest
    group.add_argument( '--maxrel',        action='store_const', const='MAXREL',                dest='MRMR_METHOD')
    group.add_argument( '--mid',           action='store_const', const='MID',                   dest='MRMR_METHOD')
    group.add_argument( '--miq',           action='store_const', const='MIQ',                   dest='MRMR_METHOD')
    parser.add_argument('--normalizemrmr', action='store_true',                                 dest='MRMR_NORMALIZE')
    parser.add_argument('--similar',                                                type=float, dest='SIMILAR')
    parser.set_defaults(
        MRMR_METHOD   ='MID',
        MRMR_NORMALIZE=False,
        SIMILAR       =0.0 # no similar features by default
        )
    return parser


def numtype(string):
    try:
        return int(string)
    except ValueError:
        try:
            return float(string)
        except ValueError:
            raise ArgumentTypeError('invalid number: "{0:s}"'.format(string))


def rfe_args(parser):
    parser.add_argument('--rfe',     action='store_true',               dest='RFE')
    parser.add_argument('--rfestep',                      type=numtype, dest='RFE_STEP')
    parser.set_defaults(
        RFE_STEP=1
        )
    return parser


def optstat_args(parser):
    group = parser.add_mutually_exclusive_group()
    #                  option           action                const                     dest
    group.add_argument('--accuracy',    action='store_const', const=Scorer.ACCURACY,    dest='OPTSTAT')
    group.add_argument('--ppv',         action='store_const', const=Scorer.PPV,         dest='OPTSTAT')
    group.add_argument('--precision',   action='store_const', const=Scorer.PPV,         dest='OPTSTAT')
    group.add_argument('--npv',         action='store_const', const=Scorer.NPV,         dest='OPTSTAT')
    group.add_argument('--sensitivity', action='store_const', const=Scorer.SENSITIVITY, dest='OPTSTAT')
    group.add_argument('--recall',      action='store_const', const=Scorer.SENSITIVITY, dest='OPTSTAT')
    group.add_argument('--specificity', action='store_const', const=Scorer.SPECIFICITY, dest='OPTSTAT')
    group.add_argument('--tnr',         action='store_const', const=Scorer.SPECIFICITY, dest='OPTSTAT')
    group.add_argument('--f1score',     action='store_const', const=Scorer.F1SCORE,     dest='OPTSTAT')
    group.add_argument('--mcc',         action='store_const', const=Scorer.MCC,         dest='OPTSTAT')
    parser.set_defaults(
        OPTSTAT=Scorer.MCC
        )
    return parser


def feature_args(parser):
    #                   option             action                type      dest
    parser.add_argument('--no-pngs',       action='store_false',           dest='PNGS')
    parser.add_argument('--no-pngs-pairs', action='store_false',           dest='PNGS_PAIRS')
    parser.add_argument('--radius',                              type=int, dest='RADIUS')
    parser.set_defaults(
        RADIUS=0
        )
    return parser


def filter_args(parser):
    #                   option      type        dest
    parser.add_argument('--maxcon', type=float, dest='MAX_CONSERVATION')
    parser.add_argument('--maxgap', type=float, dest='MAX_GAP_RATIO')
    parser.add_argument('--mincon', type=float, dest='MIN_CONSERVATION')
    parser.set_defaults(
        MAX_CONSERVATION=1.0, # 93.,
        MAX_GAP_RATIO   =0.1, # 93.,
        MIN_CONSERVATION=1.0, # 33.,
        )
    return parser


def log2ctype(string):
    try:
        begin, end, step = string.split(',')
        return (int(begin), int(end), float(step))
    except ValueError:
        raise ArgumentTypeError("'%s' is not a triple of form (int, int, float)" % string)


def svm_args(parser):
    #                   option     type          dest
    parser.add_argument('--log2c', type=log2ctype, dest='LOG2C')
    parser.set_defaults(
        LOG2C=(-5, 15, 0.25)
        )
    return parser


def cv_args(parser):
    #                   option     action               type      dest
    parser.add_argument('--cv',                         type=int, dest='CV_FOLDS')
    parser.add_argument('--loocv', action='store_true',           dest='LOOCV')
    parser.set_defaults(
        CV_FOLDS=5,
        LOOCV   =False
        )
    return parser


def simtype(string):
    if string == 'randdumbepi':
        return Simulation.DUMB
    elif string == 'randepi':
        return Simulation.EPITOPE
    elif string == 'randseq':
        return Simulation.SEQUENCE
    elif string == 'randtarget':
        return Simulation.TARGET
    else:
        raise ArgumentTypeError("'%s' is not one of %s" % (string, ', '.join(Simulation.VALUES)))


def probtype(string):
    try:
        val = float(string)
        assert 0. <= val and val <= 1.
        return val
    except:
        raise ArgumentTypeError('must be a real in range [0, 1]')


def nattype(string):
    try:
        val = int(string)
        assert 0 < val
        return val
    except:
        raise ArgumentTypeError('must be an integer in range [1..]')


def simulation_args(parser):
    #                   option             type           dest
    parser.add_argument('--sim',           type=simtype,  dest='SIM')
    parser.add_argument('--simruns',       type=nattype,  dest='SIM_RUNS')
    parser.add_argument('--simepisize',    type=nattype,  dest='SIM_EPI_SIZE')
    parser.add_argument('--simepimutrate', type=probtype, dest='SIM_EPI_MUT_RATE')
    parser.add_argument('--simepiseqnum',  type=nattype,  dest='SIM_EPI_N')
    parser.add_argument('--simepinoise',   type=probtype, dest='SIM_EPI_NOISE')
    parser.add_argument('--simepiperc',    type=probtype, dest='SIM_EPI_PERCENTILE')
    parser.set_defaults(
        SIM               =None, # can be 'randtarget' for now
        SIM_RUNS          =1,
        SIM_EPI_SIZE      =10,
        SIM_EPI_MUT_RATE  =0.01,
        SIM_EPI_N         =None, # default is to use len(seqrecords),
        SIM_EPI_NOISE     =0.08,
        SIM_EPI_PERCENTILE=0.5
        )
    return parser


def cutofftype(string):
    try:
        return float(string)
    except:
        raise ArgumentTypeError('must be a real')


class LabelTypeFactory:

    def __init__(self, data):
        self.valid_labels = sorted(data.labels)

    def __call__(self, string):
        if string not in self.valid_labels:
            msg = "'%s' is not in the list of valid labels: %s" % (
                string,
                ', '.join("'%s'" % lab for lab in self.valid_labels))
            raise ArgumentTypeError(msg)
        return string


class SubtypeTypeFactory:

    def __init__(self, data):
        self.valid_subtypes = sorted(data.subtypes, key=lambda x: x.strip().upper())

    def __call__(self, string):
        if not self.valid_subtypes:
            raise ArgumentTypeError("data source does not support subtype filtering")
        subtypes = string.split(',')
        for subtype in subtypes:
            if subtype not in self.valid_subtypes:
                msg = "'%s' not in the list of possible subtypes: %s" % (
                    subtype,
                    ', '.join("'%s'" % st.strip() for st in self.valid_subtypes))
                raise ArgumentTypeError(msg)
        return set(subtypes)


class AntibodyTypeFactory:

    def __init__(self, data):
        self.valid_antibodies = sorted(data.antibodies, key=lambda x: x.strip())

    def __call__(self, string):
        if string not in self.valid_antibodies:
            if ' ' + string not in self.valid_antibodies:
                msg = "'%s' not in the list of possible antibodies: %s" % (
                    string,
                    ', '.join("'%s'" % ab.strip() for ab in self.valid_antibodies))
                raise ArgumentTypeError(msg)
            else:
                string = ' ' + string
        return string


class FastaTypeFactory:

    def __init__(self, is_dna):
        self.is_dna = is_dna

    def __call__(self, string):
        try:
            with open(string) as h:
                source = Verifier(SeqIO.parse(h, seqfile_format(string)), DNAAlphabet)
                try:
                    seq = next(iter(source))
                    if not self.is_dna:
                        seq = translate(seq)
                except VerifyError:
                    if self.is_dna:
                        raise ArgumentTypeError("DNA encoding incompatible with protein reference")
                    source.set_alphabet(AminoAlphabet)
                    seq = next(iter(source))
            return seq
        except ArgumentTypeError:
            raise sys.exc_info()[1]
        except:
            raise ArgumentTypeError("invalid FASTA file '{0:s}'".format(string))


def PathType(string):
    from argparse import ArgumentTypeError
    if not isfile(string):
        raise ArgumentTypeError("file '%s' does not exist!" % string)
    return string


def SeedType(string):
    try:
        val = int(string)
        seed(val)
        np_seed(val)
        return val
    except ValueError:
        raise ArgumentTypeError("invalid seed '{0:s}'".format(string))


def init_args(description, args):
    from idepi import __path__ as idepi_path

    parser = ArgumentParser(description=description)

    # handle the datasource, we need to know to setup labeltype and subtype info
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--csv',    type=PathType, dest='_DATA', nargs=2, metavar=('FASTA', 'CSV'))
    group.add_argument('--sqlite', type=PathType, dest='_DATA', nargs=1, metavar='SQLITE3')
    group.set_defaults(
        _DATA=[join(idepi_path[0], 'data', 'allneuts.sqlite3')]
        )

    # handle the encoder early as well
    encoders = dict((str(enc), enc) for enc in (AminoEncoder, DNAEncoder, StanfelEncoder))
    parser.add_argument(
        '--encoding',
        type=lambda s: encoders.get(s, s),
        choices=sorted(encoders.values(), key=str),
        dest='ENCODER',
        default=AminoEncoder
        )

    # rather than removing the help and making a new parser,
    # if help options are passed defer them to the next parsing
    deferred = []
    for arg in ('-h', '--help'):
        try:
            args.remove(arg)
            deferred.append(arg)
        except ValueError:
            pass

    ns, args = parser.parse_known_args(args)

    # deferred
    args += deferred

    # setup a "subtypetype for the parser"
    is_dna = ns.ENCODER == DNAEncoder
    ns.DATA = DataSource(*ns._DATA)
    fastatype = FastaTypeFactory(is_dna)
    # labeltype = labeltypefactory(ns.DATA)
    subtype = SubtypeTypeFactory(ns.DATA)

    #                   option             action               type                dest
    parser.add_argument('--log',                                type=logtype,       dest='LOGGING')
    parser.add_argument('--label',                              type=str,           dest='LABEL')
    parser.add_argument('--filter',                             type=csvtype,       dest='FILTER')
    parser.add_argument('--clonal',        action='store_true',                     dest='CLONAL')
    parser.add_argument('--subtypes',                           type=subtype,       dest='SUBTYPES')
    parser.add_argument('--weighting',     action='store_true',                     dest='WEIGHTING')
    parser.add_argument('--refmsa',                             type=PathType,      dest='REFMSA')
    parser.add_argument('--refseq',                             type=fastatype,     dest='REFSEQ')
    parser.add_argument('--test',          action='store_true',                     dest='TEST')
    parser.add_argument('--seed',                               type=SeedType,      dest='RAND_SEED')
    parser.add_argument('-o', '--output',                       type=FileType('w'), dest='OUTPUT')

    refseq = hxb2.env.load()

    parser.set_defaults(
        LOGGING    =None,
        LABEL      ='max(IC50) > 20',
        FILTER     =[],
        CLONAL     =False,
        SUBTYPES   =set(),
        WEIGHTING  =False,
        REFMSA     =PathType(join(idepi_path[0], 'data', 'HIV1_FLT_2012_env_DNA.sto')),
        REFSEQ     =refseq if is_dna else translate(refseq),
        RAND_SEED  =42, # magic number for determinism
        PHYLOFILTER=False,
        OUTPUT     =sys.stdout
        )

    return parser, ns, args


def parse_args(parser, args, namespace=None):
    ns = parser.parse_args(args=args, namespace=namespace)
    return ns


def finalize_args(ns):
    if ns.OUTPUT != sys.stdout:
        ns.OUTPUT.close()
