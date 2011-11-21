#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
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

from operator import itemgetter
from os import close, remove, rename
from random import random
from re import sub, match
from sqlite3 import connect
from sys import stderr
from tempfile import mkstemp

import numpy as np

from Bio import SeqIO

from _abrecord import AbRecord
from _alphabet import Alphabet
from _hmmer import Hmmer
from _smldata import SmlData


__all__ = [
    'BASE_ALPH',
    'set_util_params',
    'is_testdata',
    'is_HXB2',
    'id_to_class',
    'id_to_float',
    'id_to_subtype',
    'get_noise',
    'get_valid_antibodies_from_db',
    'get_valid_subtypes_from_db',
    'base_10_to_n',
    'base_26_to_alph',
    'alph_to_base_26',
    'base_n_to_10',
    'ystoconfusionmatrix',
    'collect_AbRecords_from_db',
    'clamp',
    'sanitize_seq',
]

__HXB2_IDS = None
__IC50LT = None
__IC50GT = None

BASE_ALPH = 26


def set_util_params(hxb2_ids, ic50gt=None, ic50lt=None):
    global __HXB2_IDS, __IC50GT, __IC50LT
    __HXB2_IDS = hxb2_ids
    __IC50GT = ic50gt
    __IC50LT = ic50lt


def is_testdata(id):
    if is_HXB2(id):
        return False
    else:
        try:
            if not id.rsplit('|', 3)[3]:
                return True
            else:
                return False
        except IndexError, e:
            print >> stderr, 'ERROR: malformed ID: %s' % id
            raise e
    return False


def is_HXB2(id):
    try:
        gid = id
        gid = gid.split('|', 2)[1]
        gid = gid.split(':', 1)[0]
        if gid in __HXB2_IDS:
            return True
    except IndexError, e:
        pass
        # we don't care about this, do we?
        # print >> stderr, "ERROR: malformed ID: %s" % id
    return False


def id_to_class(id, target):
    if target not in ("lt", "gt"):
        print >> stderr, "ERROR: target must be one of \"lt\" or \"gt\"."
        raise ValueError
    try:
        ic50 = id.rsplit('|', 3)[3] # accession | subtype | ab | ic50
        if not ic50:
            return None
        else:
            ic50 = float(ic50)
    except ValueError, e:
        print >> stderr, "ERROR: cannot parse '%s' for IC50 value" % id
        raise e
    if __IC50GT is None or __IC50LT is None:
        print >> stderr, 'ERROR: call set_util_params to set IC50GT and IC50LT values for utilities'
        exit(-1)
    if target == "gt":
        c = ic50 > __IC50GT
    else:
        c = ic50 < __IC50LT
    # print '%.2f\t%d' % (ic50, c)
    return c


def id_to_subtype(id):
    try:
        subtype = id.rsplit('|', 3)[1].upper()
    except ValueError, e:
        raise ValueError('Cannot parse `%s\' for HIV subtype' % id)
    return subtype


def id_to_float(id):
    try:
        ic50 = float(id.rsplit('|', 3)[3]) # accession | subtype | ab | ic50
    except ValueError, e:
        raise ValueError('Cannot parse `%s\' for IC50 value' % id)
    return ic50


def get_noise(id):
    return id_to_float(id)


def base_10_to_n(n, N):
    val = n
    cols = {}
    pow_ = -1
    while val >= N:
        new_val = val
        mul_ = 0
        pow_ = 0
        while new_val > 0:
            new_mul = new_val / N
            if new_mul > 0:
                mul_ = new_mul
                pow_ += 1
            new_val /= N
        val -= pow(N, pow_) * mul_
        cols[pow_] = mul_
    cols[0] = val
    for i in xrange(min(cols.keys())+1, max(cols.keys())):
        if i not in cols:
            cols[i] = 0
    return cols


def base_26_to_alph(cols):
    for k in sorted(cols.keys()):
        # we might think this dangerous, but if k+1 is in cols, then it is > 1 or it has something above it
        if cols[k] <= 0 and (k+1) in cols:
            cols[k+1] -= 1
            cols[k] += 26
    if cols[max(cols.keys())] == 0:
        del cols[max(cols.keys())]
    alph = ""
    for k, v in sorted(cols.items(), key=itemgetter(0), reverse=True):
        alph += chr(ord('a') + v - 1)
    return alph


def alph_to_base_26(str):
    cols = {}
    col_idx = 0
    for i in xrange(len(str)-1, -1, -1):
        new_val = ord(str[i]) - ord('a') + 1
        cols[col_idx] = new_val
        col_idx += 1
    for i in xrange(col_idx):
        if cols[i] > 25:
            cols[i] %= 26
            if (i+1) not in cols:
                cols[i+1] = 0
            cols[i+1] += 1
    return cols


def base_n_to_10(cols, N):
    num = 0
    for k, v in cols.items():
        num += pow(N, k) * v
    return num


# very heavily based on the design of friedmanchisquare in scipy
try:
    from scipy.special import fdtrc
    def durbin(*args):

        # taken verbatim from scipy.stats._support.abut
        def _abut(source, *args):
            source = np.asarray(source)
            if len(source.shape) == 1:
                width = 1
                source = np.resize(source, [source.shape[0], width])
            else:
                width = source.shape[1]
            for addon in args:
                if len(addon.shape) == 1:
                    width = 1
                    addon = np.resize(addon, [source.shape[0], width])
                else:
                    width = source.shape[1]
                if len(addon) < len(source):
                    addon = np.resize(addon, [source.shape[0], addon.shape[1]])
                elif len(addon) > len(source):
                    source = np.resize(source, [addon.shape[0], source.shape[1]])
                source = np.concatenate((source, addon), 1)
            return source

        # also taken from scipy.stats, but ignores everything under 0.
        def _rankposdata(a):
            a = np.ravel(a)
            b = np.argsort(a)
            a = a[b]
            n = len(a)
            dupcount = 0
            oldrank = -1
            sumranks = 0
            newarray = np.zeros(n, float)
            for i in xrange(n):
                if a[i] <= 0.:
                    newarray[b[i]] = 0.
                    continue
                oldrank += 1
                sumranks += oldrank
                dupcount += 1
                if i == n-1 or a[i] != a[i+1]:
                    averrank = float(sumranks) / float(dupcount) + 1
                    for j in xrange(i-dupcount+1, i+1):
                        newarray[b[j]] = averrank
                    sumranks = 0
                    dupcount = 0
            return newarray

        b = len(args)
        if b < 3:
            raise ValueError, 'Less than 3 levels. Durbin test is not appropriate'
        k = len(args[0])
        for i in xrange(1, b):
            if len(args[i]) <> k:
                raise ValueError, 'Unequal N in durbin. Aborting.'

        data = apply(_abut,args)
        data = data.astype(float)

        A = 0.
        t = data.shape[1]
        R = np.zeros(t, float)
        rs = np.zeros(t, int)
        for i in xrange(len(data)):
            data[i] = _rankposdata(data[i])
            for j in xrange(len(data[i])):
                A += pow(data[i,j], 2.)
                R[j] += data[i,j]
                if data[i,j] > 0.:
                    rs[j] += 1

        r = np.mean(rs)
        t = float(t)
        b = float(b)
        k = float(k)
        C = b * k * pow(k + 1, 2) / 4
        T1 = (t-1) * sum(map(lambda x: pow(x, 2) - r*C, R)) / (A-C)
        T2 = (T1 / (t-1)) / ((b*k - b - T1) / (b*k - b - t + 1))

        print data
        print R
        print "r = %g, t = %g, b = %g, k = %g, C = %g, A = %g, T1 = %g" % (r, t, b, k, C, A, T1)

        return T2, fdtrc(k-1, b*k-b-t+1, T2)

    __all__ += ['durbin']

except ImportError:
    pass

def ystoconfusionmatrix(truth, preds):
    tps = truth > 0.
    tns = truth <= 0.
    pps = preds > 0.
    pns = preds <= 0.

                                                           # true pos    true neg    false pos   false neg
    tp, tn, fp, fn = map(lambda a: np.sum(np.multiply(*a)), [(tps, pps), (tns, pns), (tns, pps), (tps, pns)])

    return (tp, tn, fp, fn)

def get_valid_subtypes_from_db(dbpath):
    conn = connect(dbpath)
    curr = conn.cursor()

    curr.execute('''select distinct SUBTYPE from GENO_REPORT''')

    valid_subtypes = [r[0] for r in curr if r[0].strip() != '']

    conn.close()

    return valid_subtypes

def get_valid_antibodies_from_db(dbpath):
    conn = connect(dbpath)
    curr = conn.cursor()

    curr.execute('''select distinct ANTIBODY from NEUT''')

    valid_antibodies = [r[0] for r in curr]

    conn.close()

    return valid_antibodies

def collect_AbRecords_from_db(dbpath, antibody):
    conn = connect(dbpath)
    curr = conn.cursor()

    curr.execute('''
    select distinct S.ID as ID, S.SEQ as SEQ, G.SUBTYPE as SUBTYPE, N.AB as AB, N.IC50 as IC50 from
    (select ACCESSION_ID as ID, RAW_SEQ as SEQ from SEQUENCE group by ACCESSION_ID) as S join
    (select ACCESSION_ID as ID, SUBTYPE from GENO_REPORT group by ACCESSION_ID) as G join
    (select ACCESSION_ID as ID, ANTIBODY as AB, IC50_STRING as IC50 from NEUT where ANTIBODY = ?) as N
    on N.ID = S.ID and G.ID = S.ID order by S.ID;
    ''', (antibody,))

    # make sure the records are unique
    ab_records = []
    for row in curr:
        try:
            ab_record = AbRecord(*row)
            if ab_record.id in [abr.id for abr in ab_records]:
                ab_record.id += '-1'
            ab_records.append(ab_record)
        except ValueError:
            continue

    conn.close()

    return ab_records

def clamp(x):
    if x < 0.:
        return 0.
    if x > 1.:
        return 1.
    return x

def sanitize_seq(seq, alphabet):
    alphdict = alphabet.todict()
    assert(len(Alphabet.SPACE) > 0 and len(seq) > 0 and len(alphdict.keys()) > 0)
    try:
        seq = str(seq)
        seq = seq.upper()
        seq = sub(r'[%s]' % Alphabet.SPACE, '-', seq)
        seq = sub(r'[^%s]' % ''.join(alphdict.keys()), 'X', seq)
    except TypeError, e:
        raise RuntimeError('something is amiss with things:\n  SPACE = %s\n  seq = %s\n  alphabet = %s\n' % (Alphabet.SPACE, seq, alphdict))
    return seq
