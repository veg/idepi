#!/usr/bin/env python
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

import sys

from operator import itemgetter
from random import gauss, randint, random
from six import StringIO

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein


__all__ = ['DumbRandomSequences', 'MarkovRandomSequences']

__NAME = None
__GAP_CHARS = ('.', '_', '-', '=', 'x', 'X')
__GAP = '-'


def help():
    print('usage: %s [-N <NUMSEQ>] <ALNFILE>' % __NAME, file=sys.stderr)
    return -1


def parse_opts(argv):
    global NAME
    opts = { 'N': (int, 1), 'c': False, 'g': False }
    NAME = argv[0]
    i = 1
    while i < len(argv):
        arg = argv[i]
        if arg[0] == '-':
            v = arg[1]
            if v in opts.keys():
                if isinstance(opts[v], tuple):
                    opts[v] = (opts[v][0], opts[v][0](argv[i+1]))
                    i += 1
                else:
                    opts[v] = not opts[v]
            else:
                print('ERROR: INVALID ARGUMENT %s' % arg, file=sys.stderr)
                return None, None, help()
        else:
            return opts, argv[i:], 0
        i += 1
    return None, None, help()


# TODO: replace this idiocy with a BioPython parser...
def read_alignment(file):
    fh = open(file)

    id_seqs = {}
    seq_len = 0
    on = False

    for line in fh:
        line = line.strip()
        if line == '# STOCKHOLM 1.0':
            on = True
            continue
        elif line == '//':
            on = False
            continue
        else:
            if on:
                if line == '':
                    continue
                elif line[:2] == '#=':
                    continue
                else:
                    id, seq = line.split()
                    if id not in seq:
                        if seq != '':
                            id_seqs[id] = seq
                            if seq_len == 0:
                                seq_len = len(seq)
                        else:
                            print('ERROR: SEQUENCE IS OF LENGTH 0', file=sys.stderr)
                            return -1
                    else:
                        print('ERROR: DUPLICATE IDS IN ALIGNMENT', file=sys.stderr)
                        return -1

    fh.close()
    return seq_len, id_seqs


def seqrecords_to_alignment(seqrecords):
    buf = StringIO()
    AlignIO.write(seqrecords, buf, 'stockholm')
    buf.seek(0)
    alignment = AlignIO.read(buf, 'stockholm')
    return alignment


def DumbRandomSequences(base_seq, N=1, consensus=False, gaps=False, opts={}, idfmt='', noise=0., rate=0., alphabet=None):

    if len(opts) == 0:
        opts['N'] = (int, N)
        opts['c'] = consensus
        opts['g'] = gaps

    if alphabet is not None:
        alphabet = [a for a in alphabet if a not in ('X', '-')]

    seqrecords = []

    for n in range(opts['N'][1]):
        seq = [i for i in base_seq]
        if rate > 0.:
            for p in range(len(seq)):
                if random() < rate:
                    c = seq[p]
                    while c == seq[p]:
                        c = alphabet[randint(0, len(alphabet) - 1)]
                    seq[p] = c
        seq_str = ''.join(seq)

        # hopefully we won't need this, like, at all
        if opts['g']:
            trim_chars = (' ')
        else:
            trim_chars = (' ', __GAP)
        assert(seq_str == ''.join(seq[i] for i in range(len(seq)) if seq[i] not in trim_chars))
        seq_id = idfmt % ('sample-%i' % (n + 1), gauss(0, noise))
        seqrecords.append(SeqRecord(Seq(seq_str, generic_protein), id=seq_id, name=seq_id))

    alignment = seqrecords_to_alignment(seqrecords)

    return alignment


def MarkovRandomSequences(file, N=1, consensus=False, gaps=False, opts={}, idfmt='', noise=0., rate=0., alphabet=None):
    if rate > 0. and alphabet is None:
        print('ERROR: MarkovRandomSequences requires an alphabet if substitution rate > 0', file=sys.stderr)
        sys.exit(-1)

    if len(opts) == 0:
        opts['N'] = (int, N)
        opts['c'] = consensus
        opts['g'] = gaps

    if alphabet is not None:
        alphabet = [a for a in alphabet if a not in ('X', '-')]

    seq_len, id_seqs = read_alignment(file)

    tree = [{}] * (seq_len)

    # make the markov tree (this could be done more elegantly)
    for i in range(seq_len):
        vals = {}
        for seq in id_seqs.values():
            v = seq[i].upper()
            if v in __GAP_CHARS:
                v = __GAP
            if v not in vals:
                vals[v] = [0., {}]
            vals[v][0] += 1.
            # if we're not the last position
            if i+1 < seq_len:
                w = seq[i+1].upper()
                if w in __GAP_CHARS:
                    w = __GAP
                if w not in vals[v][1]:
                    vals[v][1][w] = 0.
                vals[v][1][w] += 1.
        # make everything into a probability distribution
        sum_v = sum(v[0] for v in vals.values())
        prev_v = 0.
        for v, _ in sorted(vals.items(), key = itemgetter(1), reverse = True):
            vals[v][0] /= sum_v
            vals[v][0] += prev_v
            prev_v = vals[v][0]
            # if we're not the last position
            if i+1 < seq_len:
                sum_w = sum(vals[v][1].values())
                prev_w = 0.
                for w, _ in sorted(vals[v][1].items(), key = itemgetter(1), reverse = True):
                    vals[v][1][w] /= sum_w
                    vals[v][1][w] += prev_w
                    prev_w = vals[v][1][w]
        tree[i] = vals

    # make the seqs
    seqrecords = []

    for n in range(opts['N'][1]):
        seq = [' '] * seq_len
        if opts['c']:
            for i in range(len(tree)):
                for w, _ in sorted(tree[i].items(), key = lambda x: x[1][0]):
                    if w == __GAP:
                        continue
                    seq[i] = w
                    break
        else:
            # get the 0th position
            r = random()
            for w, tmp_r in sorted(tree[0].items(), key = lambda x: x[1][0]):
                if tmp_r[0] > r:
                    seq[0] = w
                    break
            # get all following positions
            for i in range(len(tree)-1):
                r = random()
                for w, curr_r in sorted(tree[i][seq[i]][1].items(), key = itemgetter(1)):
                    if curr_r > r:
                        seq[i+1] = w
                        break
        # perform some base mutation on the sequences
        if rate > 0.:
            for p in range(len(seq)):
                if random() < rate:
                    c = seq[p]
                    while c == seq[p]:
                        c = alphabet[randint(0, len(alphabet) - 1)]
                    seq[p] = c
        if opts['g']:
            trim_chars = (' ')
        else:
            trim_chars = (' ', __GAP)
        seq = [seq[i] for i in range(len(seq)) if seq[i] not in trim_chars]
        seq_str = ''.join(seq)
        # seq_str = '\n'.join(['%s' % ''.join(seq[i:i+80]) for i in xrange(len(seq), 80)])
        seq_id = idfmt % ('sample-%i' % (n + 1), gauss(0, noise))
        seqrecords.append(SeqRecord(Seq(seq_str, generic_protein), id=seq_id, name=seq_id))

    alignment = seqrecords_to_alignment(seqrecords)

    return alignment


def main(argv = sys.argv):

    opts, args, ret = parse_opts(argv)
    if ret < 0:
        return ret

    seq_records = MarkovRandomSequences(args[0], opts = opts)

    # TODO: gapped output is broken
    SeqIO.write(seq_records, sys.stdout, 'fasta')

    return 0


if __name__ == '__main__':
    sys.exit(main())
