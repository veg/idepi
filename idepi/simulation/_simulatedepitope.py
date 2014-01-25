
from __future__ import division, print_function

from re import sub

from .._common import clamp, sanitize_seq


__all__ = ['SimulatedEpitope']


class SimulatedEpitope:

    def __init__(self, positions, position_names, alphabet, kernel_func=None):
        self.positions = positions
        self.names = position_names
        self.alphabet = alphabet
        # default to a uniform linear kernel
        if kernel_func is None:
            self.kernel_func = lambda x, n: clamp(1. * x / len(positions) + n)

    def __str__(self):
        return '\n'.join(sorted(self.names, key = lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x))))

    def evaluate(self, seq, noise=0., proportion=-1):
        total = 0
        # sanitize the sequence, so that it fits in our alphabet
        seq = sanitize_seq(seq, self.alphabet)
        alphdict = self.alphabet.todict()
        for k, v in self.positions.items():
            if alphdict[seq[k]] == alphdict[v]:
                total += 1
        # this thing should now produce 50/50 splits no matter what
        # if the positions are more mutated than the base rate, then 25 (resistant)
        # else 1 (susceptible)
        if proportion < 0:
            ret = self.kernel_func(total, noise)
        else:
            ret = 1. if proportion < self.kernel_func(total, noise) else 25.
        # 12.5 * pow(1. - self.kernel_func(total), 0.8) # 0.8 for a 0.1 mutation rate, 2.667 for a 50%/50% split
        # ret = abs(2. * log10(self.kernel_func(total) + 0.00018) / log10(proportion))
        return ret
