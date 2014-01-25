
from __future__ import division, print_function

from re import sub as re_sub

import numpy as np

from idepi.constants import GAPS


__all__ = [
    'BASE_ALPH',
    'base_10_to_n',
    'base_26_to_alph',
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


def get_noise(seqrecord, label='IC50'):
    from .util import seqrecord_get_values
    # just return the "mean" as noise
    return np.mean(seqrecord_get_values(seqrecord.description, label))


def sanitize_seq(seq, alphabet):
    alphdict = alphabet.todict()
    assert(len(GAPS) > 0 and len(seq) > 0 and len(alphdict) > 0)
    try:
        seq = str(seq)
        seq = seq.upper()
        seq = re_sub(r'[%s]' % GAPS, '-', seq)
        seq = re_sub(r'[^%s]' % ''.join(alphdict.keys()), 'X', seq)
    except TypeError:
        raise RuntimeError(
            'something is amiss with things:\n  GAPS = %s\n  seq = %s\n  alphabet = %s\n' % (
                GAPS, seq, alphdict)
            )
    return seq


def clamp(x):
    if x < 0.:
        return 0.
    if x > 1.:
        return 1.
    return x
