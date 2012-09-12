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

from __future__ import division, print_function

from collections import namedtuple
from os.path import exists
from re import compile as re_compile
from warnings import warn

import numpy as np

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import Gapped, generic_nucleotide, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .._common import (
    BASE_ALPH,
    base_26_to_alph,
    base_10_to_n,
    generate_alignment_from_seqrecords,
    refseq_off
)
from ..hmmer import Hmmer


__all__ = [
    'set_util_params',
    'is_refseq',
    'seqrecord_get_ic50s',
    'seqrecord_get_subtype',
    'seqrecord_set_ic50',
    'ystoconfusionmatrix',
    'alignment_identify_refidx',
    'extract_feature_weights_similar',
    'extract_feature_weights',
    'generate_alignment',
    'crude_sto_read',
    'C_range'
]

__REFSEQ_IDS = []
__IC50 = None


def set_util_params(refseq_ids=None, ic50=None):
    global __REFSEQ_IDS, __IC50
    if refseq_ids is not None:
        __REFSEQ_IDS = refseq_ids
    if ic50 is not None:
        __IC50 = ic50


def is_refseq(seqrecord):
    try:
        if seqrecord.id.strip() in __REFSEQ_IDS:
            return True
    except IndexError:
        pass
        # we don't care about this, do we?
        # print >> stderr, "ERROR: malformed ID: %s" % id
    return False


def seqrecord_get_ic50s(seqrecord):
    # cap ic50s to 25
    try:
        ic50s = [
            min(float(ic50.strip().lstrip('<>')), 25.) for ic50 in seqrecord.description.rsplit('|', 2)[2].split(',')
        ] # subtype | ab | ic50
    except ValueError:
        raise ValueError('Cannot parse `%s\' for IC50 value' % seqrecord.description)
    return ic50s


def seqrecord_get_subtype(seqrecord):
    try:
        subtype = seqrecord.description.rsplit('|', 2)[0].upper()
    except ValueError:
        raise ValueError('Cannot parse `%s\' for HIV subtype' % seqrecord.description)
    return subtype


def seqrecord_set_ic50(seqrecord, ic50):
    vals = seqrecord.description.rsplit('|', 2)
    while len(vals) < 3:
        vals.append('')
    vals[2] = str(ic50)
    seqrecord.description = '|'.join(vals)
    return seqrecord


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
            for i in range(n):
                if a[i] <= 0.:
                    newarray[b[i]] = 0.
                    continue
                oldrank += 1
                sumranks += oldrank
                dupcount += 1
                if i == n-1 or a[i] != a[i+1]:
                    averrank = float(sumranks) / float(dupcount) + 1
                    for j in range(i-dupcount+1, i+1):
                        newarray[b[j]] = averrank
                    sumranks = 0
                    dupcount = 0
            return newarray

        b = len(args)
        if b < 3:
            raise ValueError('Less than 3 levels. Durbin test is not appropriate')
        k = len(args[0])
        for i in range(1, b):
            if len(args[i]) != k:
                raise ValueError('Unequal N in durbin. Aborting.')

        data = _abut(*args)
        data = data.astype(float)

        A = 0.
        t = data.shape[1]
        R = np.zeros(t, float)
        rs = np.zeros(t, int)
        for i in range(len(data)):
            data[i] = _rankposdata(data[i])
            for j in range(len(data[i])):
                A += pow(data[i,j], 2.)
                R[j] += data[i,j]
                if data[i,j] > 0.:
                    rs[j] += 1

        r = np.mean(rs)
        t = float(t)
        b = float(b)
        k = float(k)
        C = b * k * pow(k + 1, 2) / 4
        T1 = (t-1) * sum([pow(x, 2) - r*C for x in R]) / (A-C)
        T2 = (T1 / (t-1)) / ((b*k - b - T1) / (b*k - b - t + 1))

        print(data)
        print(R)
        print("r = %g, t = %g, b = %g, k = %g, C = %g, A = %g, T1 = %g" % (r, t, b, k, C, A, T1))

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
    tp, tn, fp, fn = (np.sum(np.multiply(a, b)) for a, b in ((tps, pps), (tns, pns), (tns, pps), (tps, pns)))
    return np.array([[tp, fn], [fp, tn]], dtype=int)


def alignment_identify_refidx(alignment, ref_id_func=None):
    refidx = None
    if ref_id_func is not None:
        for i, seq in enumerate(alignment):
            if ref_id_func(seq):
                refidx = i
                break
        if refidx is None:
            raise RuntimeError('ref_id_func provided but no reference found')
    return refidx


def extract_feature_weights_similar(instance, similar=True):
    ret = {
        'features': instance.features(),
        'weights':  instance.classifier.weights()
    }
    if similar:
        ret['similar'] = instance.selector.related()
    return ret


def extract_feature_weights(instance):
    return extract_feature_weights_similar(instance, False)


def generate_alignment(seqrecords, my_basename, ref_id_func, opts):
    from ..simulation import Simulation

    sto_filename = my_basename + '.sto'
    hmm_filename = my_basename + '.hmm'

    if hasattr(opts, 'SIM') and opts.SIM == Simulation.DUMB:
        # we're assuming pre-aligned because they're all generated from the same refseq
        with open(hmm_filename, 'w') as fh:
            SeqIO.write(seqrecords, fh, 'stockholm')
    elif not exists(sto_filename):
        generate_alignment_from_seqrecords(
            seqrecords,
            my_basename,
            opts
        )

    if not exists(hmm_filename):
        hmmer = Hmmer(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)
        hmmer.build(hmm_filename, sto_filename)

    with open(sto_filename) as fh:
        msa = AlignIO.read(fh, 'stockholm')

    # we store the antibody information in the description, so grab it
    return msa, {} # crude_sto_read(sto_filename, ref_id_func, opts.DNA, description=True)


def crude_sto_read(filename, ref_id_func=None, dna=False, description=False):
    Fake = namedtuple('FakeSeqRecord', ['id'])
    alph = Gapped(generic_nucleotide if dna else generic_protein)
    refseq = None
    msa = MultipleSeqAlignment([], alphabet=alph)
    with open(filename) as fh:
        notrel = re_compile(r'^(?://|#(?!=GS))')
        isdesc = re_compile(r'^#=GS')
        trim = re_compile(r'[^-A-Z]+')
        descs = {}
        for line in fh:
            line = line.strip()
            if line == '' or notrel.match(line):
                continue
            elif isdesc.match(line):
                if description:
                    elems = line.split(None, 3)
                    if len(elems) > 3 and elems[2] == 'DE':
                        acc, desc = elems[1], elems[3]
                        if acc in descs:
                            warn("duplicate sequence name '%s' detected! The stockholm specification doesn't allow this!" % acc)
                        descs[acc] = desc
                else:
                     continue
            else:
                try:
                    acc, seq = line.split(None, 1)
                except ValueError:
                    warn("skipping line '%s', doesn't seem to contain a sequence" % line)
                    continue
                if ref_id_func is not None and ref_id_func(Fake(acc)):
                    refseq = seq
                try:
                    seq = trim.sub('', seq)
                    desc = descs[acc] if acc in descs else acc
                    msa.append(SeqRecord(Seq(seq), id=acc, description=desc))
                except ValueError:
                    warn("skipping sequence '%s', it doesn't match the length of the MSA (%d vs %d)" % (acc, len(seq), msa.get_alignment_length()))

    if ref_id_func is not None and refseq is None:
        raise RuntimeError('Unable to find the reference sequence to compute an offset!')

    offs = None if ref_id_func is None else refseq_off(refseq)
    return msa, offs


def site_labels(refseq, refseq_offs={}):
    if isinstance(refseq, str):
        ref = refseq
    elif isinstance(refseq, SeqRecord):
        ref = str(refseq.seq)
    colnames = []
    colnum = 0
    insert = 0
    for i, p in enumerate(ref):
        if i in refseq_offs:
            colnum += refseq_offs[i]
        if ref[i] not in '._-':
            colnum += 1
            insert = 0
        else:
            insert += 1
        colname = '%s%d%s' % (p if insert == 0 else '', colnum, base_26_to_alph(base_10_to_n(insert, BASE_ALPH)))
        colnames.append(colname)
    return colnames

def C_range(begin, end, step):
    recip = 1
    if isinstance(step, float):
        recip = 1. / step
        begin, end = int(recip * begin), int(recip * end)
        step = 1
    return [pow(2., float(c) / recip) for c in range(begin, end + 1, step)]
