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
from _seqtable import SeqTable
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
    'binarize_row',
    'sanitize_seq',
    'generate_relevant_data',
    'generate_feature_names',
    'compute_relevant_features_and_names',
    'binarize',
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
        print >> stderr, "ERROR: malformed ID: %s" % id
        raise e
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
    ab_records = list()
    for row in curr:
        ab_record = AbRecord(*row)
        if ab_record.id in [abr.id for abr in ab_records]:
            ab_record.id += '-1'
        ab_records.append(ab_record)

    conn.close()

    return ab_records

def clamp(x):
    if x < 0.:
        return 0.
    if x > 1.:
        return 1.
    return x

def binarize_row(row, alphabet):
    alphdict = alphabet.todict()
    alphabet_len = len(set(alphdict.values()))
    row = sanitize_seq(str(row), alphabet)
    ret = []
    for p in row:
        ret.extend([1 if i == alphdict[p] else 0 for i in xrange(alphabet_len)])
    return ret

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

def generate_relevant_data(feature_names, seq_table, alphabet, feature_idxs, target, subtypes=[], simulation=None):
    data = SmlData(feature_names)

    neg = 0
    pos = 0

    for row in seq_table.rows():
        # everything before the continue takes care of the | and : separators
        if is_HXB2(row.id):
            continue
        if len(subtypes) > 0:
            rowtype = id_to_subtype(row.id)
            if rowtype == '' or len([i for i in rowtype if i in subtypes]) <= 0:
                raise ValueError('Unwanted subtypes should already be masked out: %s' % row.id)
        expanded_row = binarize_row(row.seq, alphabet)
        feats = dict([f for f in zip(range(0, len(feature_idxs)), [expanded_row[i] for i in feature_idxs]) if f[1] != 0])
        if simulation is not None:
            class_ = random()
            class_ = 1 if class_ < simulation.proportion else 0
        else:
            class_ = id_to_class(row.id, target)
        data.add(class_, feats)
        if class_:
            pos += 1
        else:
            neg += 1

    # silence this during testing
    # print >> sys.stdout, 'Ratio + to - : %d to %d' % (pos, neg)

    return data

def generate_feature_names(seq_table, alph):
    alphabet, alphabet_names = alph.todict(), alph.names()

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
        if p not in Alphabet.SPACE:
            c += 1
            ins = 0
        else:
            ins += 1
        for v in set(alphabet.values()):
            insert = base_26_to_alph(base_10_to_n(ins, BASE_ALPH))
            names.append('%s%d%s%s' % ('' if insert != '' else p.upper(), c, insert, alphabet_names[v]))

    return names

def compute_relevant_features_and_names(seq_table, alph, names, limits={ 'max': 1.0, 'min': 1.0, 'gap': 0.1 }, filter_list=[]):
    if type(limits) != dict:
        raise ValueError('limits must be of type `dict\'')

    for lim in ('max', 'min', 'gap'):
        if lim not in limits:
            raise ValueError('limits must contain `%s\'' % lim)

    if not seq_table.loaded:
        seq_table.fill_columns()

    columns = range(0, seq_table.num_columns)
    alphabet, alphabet_names = alph.todict(), alph.names()
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
        if max_count > limits['max'] or \
           min_count > limits['min'] or \
           gap_count > limits['gap']:
            delete_cols.extend(range(j, j+alphabet_len))
            continue
        max_counts.append((i, limits['max'] - max_count))
        min_counts.append((i, limits['min'] - min_count))
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

    if len(filter_list) != 0:
        delete_cols = []
        for i in xrange(0, len(columns)):
            m = match(r'^([A-Z]\d+)', trimmed_names[i])
            if m:
                if m.group(1) not in filter_list:
                    delete_cols.append(i)
        for i in sorted(delete_cols, reverse=True):
            del columns[i]
            del trimmed_names[i]

    return (columns, trimmed_names)


def binarize(x, colnames=None, dox=True):
    assert(x.dtype == int or x.dtype == bool)
    if colnames is None and not dox:
        raise RuntimeError('binarize() does no work under these parameters!')
    nrow, ncol = x.shape
    donames = colnames is not None
    if donames:
        assert(ncol == len(colnames))
    coldim = np.zeros((ncol,), dtype=int)
    for j in xrange(ncol):
        m = max(x[:, j]) + 1
        coldim[j] = m if m > 2 else 1 if m > 1 else 0
    newcol = np.sum(coldim)
    if dox:
        newx = np.zeros((nrow, newcol), dtype=bool)
    if donames:
        newcolnames = []
    i = 0
    for j in xrange(ncol):
        if coldim[j] == 1:
            continue
        elif coldim[j] == 2:
            if dox:
                newx[:, i] = x[:, j] > 0 # binary 1 stays binary 1
            if donames:
                newcolnames.append(colnames[j])
            i += 1
        else:
            for k in xrange(coldim[j]):
                if dox:
                    newx[:, i] = x[:, j] == k
                if donames:
                    newcolnames.append(colnames[j] + '-' + base_26_to_alph(base_10_to_n(k, 26)).upper())
                i += 1
    if dox and donames:
        return newx, newcolnames
    elif dox:
        return newx
    elif donames:
        return newcolnames
    else:
        assert(0) # dead block, I think
