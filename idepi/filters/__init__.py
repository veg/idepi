
from __future__ import division, print_function

from functools import partial

from BioExt.collections import Counter

from idepi.constants import GAPS


__all__ = ['naive_filter', 'no_filter']


def __max_min_sum(iterable):
    # M == max and m == min
    M = m = next(iterable, None)
    s = 0
    if M is not None:
        s += M
        for v in iterable:
            if v > M:
                M = v
            if v < m:
                m = v
            s += v
    return M, m, s


def __naive_filter(max_conservation, min_conservation, max_gap_ratio, column):
    counts = Counter(l.upper() for l in column)
    gap_ratio = sum(counts[g] for g in GAPS) / sum(counts.values())
    maxc, minc, total = __max_min_sum(v for k, v in counts.items() if k not in GAPS)
    if (gap_ratio > max_gap_ratio or
            maxc / total > max_conservation or
            minc / total > min_conservation):
        return []
    else:
        return sorted(counts.keys())


def naive_filter(max_conservation, min_conservation, max_gap_ratio):
    return partial(__naive_filter, max_conservation, min_conservation, max_gap_ratio)


def null_filter(column):
    return sorted(set(l.upper() for l in column))
