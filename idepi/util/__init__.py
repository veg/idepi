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

from json import dumps as json_dumps, loads as json_loads
from logging import getLogger
from math import copysign
from os import close, remove
from os.path import exists, splitext
from re import compile as re_compile, I as re_I
from shutil import copyfile
from tempfile import mkstemp

import numpy as np

from Bio import AlignIO, SeqIO

from BioExt.misc import translate

from idepi.constants import AminoAlphabet, DNAAlphabet
from idepi.encoder import DNAEncoder
from idepi.hmmer import HMMER
from idepi.labeledmsa import LabeledMSA
from idepi.logging import IDEPI_LOGGER
from idepi.verifier import VerifyError, Verifier


__all__ = [
    'seqfile_format',
    'set_util_params',
    'is_refseq',
    'seqrecord_get_values',
    'seqrecord_get_subtype',
    'seqrecord_set_values',
    'ystoconfusionmatrix',
    'reference_index',
    'extract_feature_weights_similar',
    'extract_feature_weights',
    'generate_alignment',
    'C_range',
    'load_stockholm',
    'coefs_ranks'
]

__REFSEQ_IDS = []


def set_util_params(refseq_ids=None):
    global __REFSEQ_IDS
    if refseq_ids is not None:
        if isinstance(refseq_ids, (list, tuple)):
            __REFSEQ_IDS = refseq_ids
        else:
            __REFSEQ_IDS = [refseq_ids]


def is_refseq(seqrecord):
    try:
        if seqrecord.id.strip() in __REFSEQ_IDS:
            return True
    except IndexError:
        pass
        # we don't care about this, do we?
        # print >> stderr, "ERROR: malformed ID: %s" % id
    return False


def seqrecord_get_values(seqrecord, label='IC50'):
    # cap ic50s to 25
    try:
        values = json_loads(seqrecord.description)['values'][label]
    except ValueError:
        raise ValueError("Cannot parse `{0}' for {1} value".format(
            seqrecord.description,
            label
            ))
    except KeyError:
        return None
    return values


def seqrecord_get_subtype(seqrecord):
    try:
        subtype = json_loads(seqrecord.description)['subtype']
    except ValueError:
        raise ValueError("Cannot parse `%s' for HIV subtype".format(
            seqrecord.description
            ))
    return subtype


def seqrecord_set_values(seqrecord, label, values):
    desc = json_loads(seqrecord.description)
    desc['values'][label] = values
    seqrecord.description = json_dumps(desc)
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
                if i == n - 1 or a[i] != a[i + 1]:
                    averrank = float(sumranks) / float(dupcount) + 1
                    for j in range(i - dupcount + 1, i + 1):
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
                A += pow(data[i, j], 2.)
                R[j] += data[i, j]
                if data[i, j] > 0.:
                    rs[j] += 1

        r = np.mean(rs)
        t = float(t)
        b = float(b)
        k = float(k)
        C = b * k * pow(k + 1, 2) / 4
        T1 = (t - 1) * sum([pow(x, 2) - r * C for x in R]) / (A - C)
        T2 = (T1 / (t - 1)) / ((b * k - b - T1) / (b * k - b - t + 1))

        print(data)
        print(R)
        print("r = %g, t = %g, b = %g, k = %g, C = %g, A = %g, T1 = %g" % (r, t, b, k, C, A, T1))

        return T2, fdtrc(k - 1, b * k - b - t + 1, T2)

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


def reference_index(alignment, ref_id_func=None):
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
        'weights': instance.classifier.weights()
    }
    if similar:
        ret['similar'] = instance.selector.related()
    return ret


def extract_feature_weights(instance):
    return extract_feature_weights_similar(instance, False)


def seqfile_format(filename):
    return 'stockholm' if splitext(filename)[1].find('sto') == 1 else 'fasta'


def generate_hmm_(opts):
    fd, tmphmm = mkstemp(); close(fd)
    fd, tmpaln = mkstemp(); close(fd)

    is_dna = opts.ENCODER == DNAEncoder

    try:
        with open(opts.REFMSA) as msa_fh:
            with open(tmpaln, 'w') as aln_fh:
                msa_fmt = seqfile_format(opts.REFMSA)
                source = Verifier(SeqIO.parse(msa_fh, msa_fmt), DNAAlphabet)
                try:
                    SeqIO.write(
                        (record if is_dna else translate(record) for record in source),
                        aln_fh,
                        'stockholm')
                except VerifyError:
                    if is_dna:
                        raise RuntimeError("DNA encoding incompatible with protein reference MSA")
                    source.set_alphabet(AminoAlphabet)
                    aln_fh.seek(0)
                    SeqIO.write(
                        source,
                        aln_fh,
                        'stockholm')

        hmmer = HMMER(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)
        hmmer.build(
            tmphmm,
            tmpaln,
            alphabet=HMMER.DNA if is_dna else HMMER.AMINO
            )
    finally:
        if exists(tmpaln):
            remove(tmpaln)

    return tmphmm


def generate_alignment_(seqrecords, hmmfile, opts, refseq=None):
    fd, tmpseq = mkstemp(); close(fd)
    fd, tmpaln = mkstemp(); close(fd)
    finished = False

    is_dna = opts.ENCODER == DNAEncoder
    log = getLogger(IDEPI_LOGGER)

    try:
        # get the FASTA format file so we can HMMER it
        with open(tmpseq, 'w') as seq_fh:

            def records():
                if refseq:
                    yield HMMER.valid(refseq, is_dna=is_dna)
                for record in seqrecords:
                    if not is_dna and record.seq.alphabet == DNAAlphabet:
                        record = translate(record)
                    yield HMMER.valid(record, is_dna=is_dna)

            SeqIO.write(records(), seq_fh, 'fasta')

        log.debug('aligning sequences')

        hmmer = HMMER(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)
        hmmer.align(
            hmmfile,
            tmpseq,
            output=tmpaln,
            alphabet=HMMER.DNA if is_dna else HMMER.AMINO,
            outformat=HMMER.PFAM
        )

        # rename the final alignment to its destination
        finished = True
    finally:
        # cleanup these files
        if exists(tmpseq):
            remove(tmpseq)

    if not finished:
        raise RuntimeError("failed to generate alignment")

    return tmpaln


def generate_alignment(seqrecords, sto_filename, ref_id_func, opts, load=True):
    from ..simulation import Simulation

    log = getLogger(IDEPI_LOGGER)
    hmm = None

    if hasattr(opts, 'SIM') and opts.SIM == Simulation.DUMB:
        # we're assuming pre-aligned because they're all generated from the same refseq
        with open(sto_filename, 'w') as fh:
            SeqIO.write(seqrecords, fh, 'stockholm')
    else:
        try:
            tmphmm = generate_hmm_(opts)
            tmpaln = generate_alignment_(seqrecords, tmphmm, opts, refseq=opts.REFSEQ)
            copyfile(tmpaln, sto_filename)
            log.debug('finished alignment, output moved to {0:s}'.format(sto_filename))
            with open(tmphmm, 'rb') as hmm_fh:
                hmm = hmm_fh.read()
        finally:
            if exists(tmphmm):
                remove(tmphmm)
            if exists(tmpaln):
                remove(tmpaln)

    if load:
        with open(sto_filename) as fh:
            msa = AlignIO.read(fh, 'stockholm')
        refidx = reference_index(msa, ref_id_func)
        msa = LabeledMSA.from_msa_with_ref(msa, refidx)
        ranges = stockholm_rf_ranges(sto_filename)
        return trim_msa_to_ranges(msa, ranges), hmm

    return None, hmm


def C_range(begin, end, step):
    recip = 1
    if isinstance(step, float):
        recip = 1. / step
        begin, end = int(recip * begin), int(recip * end)
        step = 1
    return [pow(2., float(c) / recip) for c in range(begin, end + 1, step)]


def stockholm_rf_ranges(filename):
    hdr = re_compile('^#=GC\s+RF\s+(.+)$', re_I)
    keep = []
    with open(filename) as h:
        for line in h:
            m = hdr.match(line)
            if m:
                keep = [l not in '.-~' for l in m.group(1)]
                break
    ranges = []
    lwr = 0
    val = keep[lwr]
    for i, v in enumerate(keep):
        if v != val:
            # transition T->F, so append range
            if val:
                ranges.append((lwr, i))
            # transition F->T, so update lower bound
            else:
                lwr = i
            # update val
            val = v
    # if we have a terminal val, append final range
    if val:
        ranges.append((lwr, len(keep)))
    return ranges


def trim_msa_to_ranges(msa, ranges):
    ranges_ = iter(ranges)
    lwr, upr = next(ranges_)
    msa2 = msa[:, lwr:upr]
    for lwr, upr in ranges_:
        msa2 += msa[:, lwr:upr]
    return msa2


def load_stockholm(filename, trim=False):
    with open(filename) as h:
        msa = AlignIO.read(h, 'stockholm')
    if trim:
        ranges = stockholm_rf_ranges(filename)
        msa = trim_msa_to_ranges(msa, ranges)
    return msa


def coefs_ranks(ranking_, support_, coef_):
    coefs = {}
    ranks = {}
    col = 0
    for i, (rank, selected) in enumerate(zip(ranking_, support_)):
        if selected:
            coefs[i] = int(copysign(1, coef_[0, col]))
            ranks[i] = int(rank)
            col += 1
    return coefs, ranks
