
from __future__ import division, print_function

from logging import getLogger
from os import close, remove
from re import sub as re_sub
from shutil import copyfile
from tempfile import mkstemp

from Bio import SeqIO
from BioExt.misc import translate

from .alphabet import Alphabet
from .hmmer import HMMER
from .logging import IDEPI_LOGGER

import numpy as np


__all__ = [
    'BASE_ALPH',
    'base_10_to_n',
    'base_26_to_alph',
    'generate_alignment_from_seqrecords',
    'get_noise',
    'sanitize_seq',
    'clamp'
]

BASE_ALPH = 26


def base_10_to_n(n, N):
    if n < 0:
        sign = -1
    elif n == 0:
        return [0]
    else:
        sign = 1
    n *= sign
    digits = []
    while n:
        digits.append(n % N)
        n //= N
    return digits


def base_26_to_alph(cols):
    for i, v in enumerate(cols):
        if v <= 0 and (i + 1) < len(cols):
            cols[i + 1] -= 1
            cols[i] += 26
    if cols[-1] == 0:
        cols.pop()
    alph = ''
    for v in reversed(cols):
        alph += chr(ord('a') + v - 1)
    return alph


# def alph_to_base_26(str):
#     cols = {}
#     col_idx = 0
#     for i in range(len(str)-1, -1, -1):
#         new_val = ord(str[i]) - ord('a') + 1
#         cols[col_idx] = new_val
#         col_idx += 1
#     for i in range(col_idx):
#         if cols[i] > 25:
#             cols[i] %= 26
#             if (i+1) not in cols:
#                 cols[i+1] = 0
#             cols[i+1] += 1
#     return cols


# def base_n_to_10(cols, N):
#     num = 0
#     for k, v in cols.items():
#         num += pow(N, k) * v
#     return num


def generate_alignment_from_seqrecords(seq_records, my_basename, opts):
    fd, ab_fasta_filename = mkstemp(); close(fd)
    fd, hmm_filename = mkstemp(); close(fd)
    fd, sto_filename = mkstemp(); close(fd)
    finished = False

    log = getLogger(IDEPI_LOGGER)

    try:
        # get the FASTA format file so we can HMMER it
        with open(ab_fasta_filename, 'w') as fafh:

            def records():
                yield HMMER.valid(opts.REFSEQ)
                for record in seq_records:
                    yield HMMER.valid(record)

            SeqIO.write(records(), fafh, 'fasta')

        with open(opts.REFMSA) as msa_fh:
            with open(sto_filename, 'w') as sto_fh:
                SeqIO.write(
                    (r if opts.DNA else translate(r) for r in SeqIO.parse(msa_fh, 'stockholm')),
                    sto_fh,
                    'stockholm'
                    )

        hmmer = HMMER(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)
        numseqs = len(seq_records)

        log.debug('aligning {0:d} sequences'.format(numseqs))
        hmmer.build(hmm_filename, sto_filename, alphabet=HMMER.DNA if opts.DNA else HMMER.AMINO)
        hmmer.align(
            hmm_filename,
            ab_fasta_filename,
            output=sto_filename,
            alphabet=HMMER.DNA if opts.DNA else HMMER.AMINO,
            outformat=HMMER.PFAM
        )

        # rename the final alignment to its destination
        finished = True

    finally:
        # cleanup these files
        if finished:
            copyfile(sto_filename, my_basename + '.sto')
            log.debug('finished alignment, output moved to %s.sto' % my_basename)
        remove(sto_filename)
        remove(ab_fasta_filename)
        remove(hmm_filename)


def get_noise(seqrecord, label='IC50'):
    from .util import seqrecord_get_values
    # just return the "mean" as noise
    return np.mean(seqrecord_get_values(seqrecord.description, label))


def sanitize_seq(seq, alphabet):
    alphdict = alphabet.todict()
    assert(len(Alphabet.GAPS) > 0 and len(seq) > 0 and len(alphdict) > 0)
    try:
        seq = str(seq)
        seq = seq.upper()
        seq = re_sub(r'[%s]' % Alphabet.GAPS, '-', seq)
        seq = re_sub(r'[^%s]' % ''.join(alphdict.keys()), 'X', seq)
    except TypeError:
        raise RuntimeError(
            'something is amiss with things:\n  GAPS = %s\n  seq = %s\n  alphabet = %s\n' % (
                Alphabet.GAPS, seq, alphdict)
            )
    return seq


def clamp(x):
    if x < 0.:
        return 0.
    if x > 1.:
        return 1.
    return x
