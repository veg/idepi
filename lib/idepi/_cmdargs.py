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
#

from __future__ import division, print_function

import logging

from argparse import ArgumentParser, ArgumentTypeError
from os.path import abspath, join, split
from warnings import warn

from BioExt import hxb2

from mrmr import MRMR_LOGGER, DiscreteMrmr as MRMR
from pyxval import DiscretePerfStats as DPS

from ._alphabet import Alphabet as ALPH
from ._data import input_data as datatype
from ._simulation import Simulation


def csvtype(string):
    return string.split(',')

def logtype(string):
    loggers = set(v.lower() for v in string.split(','))
    all = set(['idepi', 'mrmr', 'pyxval', 'fakemp'])
    if 'all' in loggers:
        loggers |= all
    diff = all - loggers
    if len(diff):
        raise ValueError("unknown loggers requested: %s" % ', '.join("'%s'" for l in diff))
    # set the loggers
    if 'idepi' in loggers:
        from ._logging import IDEPI_LOGGER
        logging.getLogger(IDEPI_LOGGER).setLevel(logging.DEBUG)
    if 'mrmr' in loggers:
        logging.getLogger(MRMR_LOGGER).setLevel(logging.DEBUG)
    if 'pyxval' in loggers:
        from pyxval import PYXVAL_LOGGER
        logging.getLogger(PYXVAL_LOGGER).setLevel(logging.DEBUG)
    if 'fakemp' in loggers:
        try:
            from fakemp import FAKEMP_LOGGER
            logging.getLogger(FAKEMP_LOGGER).setLevel(logging.DEBUG)
        except ImportError:
            warn('fakemp logger not found')
    return loggers

def hmmer_args(parser):
    #                   option        type      dest
    parser.add_argument('--hmmalign', type=str, dest='HMMER_ALIGN_BIN')
    parser.add_argument('--hmmbuild', type=str, dest='HMMER_BUILD_BIN')
    parser.add_argument('--hmmiter',  type=int, dest='HMMER_ITER')
    parser.set_defaults(
        HMMER_ALIGN_BIN='hmmalign',
        HMMER_BUILD_BIN='hmmbuild',
        HMMER_ITER     =8
    )
    return parser

def featsel_args(parser):
    #                   option             action                type      dest
    parser.add_argument('--numfeats',                            type=int, dest='NUM_FEATURES')
    parser.add_argument('--forward',       action='store_true',            dest='FORWARD_SELECT')
    parser.set_defaults(
        NUM_FEATURES  =10,
        FORWARD_SELECT=False
    )
    return parser

def mrmr_args(parser):
    group = parser.add_mutually_exclusive_group()
    #                   option             action                const              dest
    group.add_argument( '--maxrel',        action='store_const', const=MRMR.MAXREL, dest='MRMR_METHOD')
    group.add_argument( '--mid',           action='store_const', const=MRMR.MID,    dest='MRMR_METHOD')
    group.add_argument( '--miq',           action='store_const', const=MRMR.MIQ,    dest='MRMR_METHOD')
    parser.add_argument('--normalizemrmr', action='store_true',                     dest='MRMR_NORMALIZE')
    parser.set_defaults(
        MRMR_METHOD   =MRMR.MID,
        MRMR_NORMALIZE=False,
        MAXREL        =False
    )
    return parser

def optstat_args(parser):
    group = parser.add_mutually_exclusive_group()
    #                  option           action                const                  dest
    group.add_argument('--accuracy',    action='store_const', const=DPS.ACCURACY,    dest='OPTSTAT')
    group.add_argument('--ppv',         action='store_const', const=DPS.PPV,         dest='OPTSTAT')
    group.add_argument('--precision',   action='store_const', const=DPS.PPV,         dest='OPTSTAT')
    group.add_argument('--npv',         action='store_const', const=DPS.NPV,         dest='OPTSTAT')
    group.add_argument('--sensitivity', action='store_const', const=DPS.SENSITIVITY, dest='OPTSTAT')
    group.add_argument('--recall',      action='store_const', const=DPS.SENSITIVITY, dest='OPTSTAT')
    group.add_argument('--specificity', action='store_const', const=DPS.SPECIFICITY, dest='OPTSTAT')
    group.add_argument('--tnr',         action='store_const', const=DPS.SPECIFICITY, dest='OPTSTAT')
    group.add_argument('--fscore',      action='store_const', const=DPS.FSCORE,      dest='OPTSTAT')
    parser.set_defaults(
        OPTSTAT=DPS.MINSTAT
    )
    return parser

def encoding_args(parser):
    group = parser.add_mutually_exclusive_group()
    #                  option       action                const               dest
    group.add_argument('--amino',   action='store_const', const=ALPH.AMINO,   dest='ALPHABET')
    group.add_argument('--dna',     action='store_const', const=ALPH.DNA,     dest='ALPHABET')
    group.add_argument('--stanfel', action='store_const', const=ALPH.STANFEL, dest='ALPHABET')
    parser.set_defaults(
        ALPHABET=ALPH.AMINO
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
        MIN_CONSERVATION=1.0  # 33.
    )
    return parser

def log2ctype(string):
    try:
        start, end, step = string.split(',')
        return (int(start), int(end), float(step))
    except ValueError:
        ArgumentTypeError("'%s' is not a triple of form (int, int, float)" % string)

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
        ArgumentTypeError("'%s' is not one of %s" % (string, ', '.join(Simulation.VALUES)))

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

def ic50type(string):
    try:
        val = float(string)
        assert 0 < val and val < 25
        return val
    except:
        raise ArgumentTypeError('must be a real in range (0, 25)')

def init_args(description):

    parser = ArgumentParser(description)

    #                   option             action='store'       type           dest
    parser.add_argument('--log',                                type=logtype,  dest='LOGGING')
    parser.add_argument('--filter',                             type=csvtype,  dest='FILTER')
    parser.add_argument('--clonal',        action='store_true',                dest='CLONAL')
    parser.add_argument('--subtypes',                           type=csvtype,  dest='SUBTYPES')
    parser.add_argument('--weighting',     action='store_true',                dest='WEIGHTING')
    parser.add_argument('--ic50',                               type=ic50type, dest='IC50')
    parser.add_argument('--data',                               type=datatype, dest='DATA')
    parser.add_argument('--refseq',                             type=str,      dest='REFSEQ')
    parser.add_argument('--ids',                                type=csvtype,  dest='REFSEQ_IDS')
    parser.add_argument('--test',          action='store_true',                dest='TEST')
    parser.add_argument('--seed',                               type=int,      dest='RAND_SEED')
    parser.add_argument('--phylofilt',     action='store_true',                dest='PHYLOFILTER')
    parser.add_argument('-o', '--output',                       type=str,      dest='OUTPUT')

    basedir = split(abspath(__file__))[0]
    refseq = hxb2.env.load()

    parser.set_defaults(
        LOGGING    =None,
        FILTER     =[],
        CLONAL     =False,
        SUBTYPES   =[],
        WEIGHTING  =False,
        IC50       =20.,
        DATA       =datatype(join(basedir, 'data', 'allneuts.sqlite3')),
        REFSEQ     =refseq,
        REFSEQ_IDS =[refseq.id],
        RAND_SEED  =42, # magic number for determinism
        PHYLOFILTER=False,
        OUTPUT     =None
    )

    return parser

def finalize_args(args):
    setattr(args, 'DNA', args.ALPHABET == Alphabet.DNA)
    return args
