
from collections import namedtuple

from BioExt.collections import Counter

from idepi.labeledmsa import LabeledMSA


PosData = namedtuple('PosData', [
    'counts',
    'gapcount',
    'label',
    'maxcount',
    'mincount',
    'pos',
    'total'
])


def posstream(alignment, alphabet):
    ncol = alignment.get_alignment_length()

    if isinstance(alignment, LabeledMSA):
        labels = zip(alignment.positions, alignment.labels)
    else:
        raise TypeError("invalid msa type: {0:s}".format(type(alignment)))

    # save the gap coord for later use
    gap = alphabet('-')

    counters = [Counter() for _ in range(ncol)]

    for i, seq in enumerate(alignment):
        for j, char in enumerate(seq):
            # update() requires an iterable, make it a 1-tuple
            char_ = (alphabet(char.upper()),)
            counters[j].update(char_)

    for counts in counters:
        pos, label = next(labels)
        gapc = 0 if gap not in counts else counts[gap]
        maxc, minc, total = maxminsum(v for k, v in counts.items() if k != gap)
        # if everything is a gap (ie, total == 0), skip it
        if total == 0:
            yield None
        else:
            yield PosData(counts, gapc, label, maxc, minc, pos, total)


def maxminsum(iterable):
    # M == max and m == min
    M = m = next(iterable, None)
    sum = 0
    if M is not None:
        sum += M
        for v in iterable:
            if v > M:
                M = v
            if v < m:
                m = v
            sum += v
    return M, m, sum
