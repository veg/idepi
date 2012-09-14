
from collections import namedtuple

from BioExt import Counter

from ..alphabet import Alphabet


PosData = namedtuple('PosData', [
    'counts',
    'gapcount',
    'insert',
    'label',
    'maxcount',
    'mincount',
    'pos',
    'total'
])


def posstream(alignment, alphabet, refidx):
    nrow = len(alignment)

    if refidx < 0 or refidx >= nrow:
        raise ValueError('refidx must be an index into alignment')

    labels = __poslabels(alignment, alphabet, refidx)

    # save the gap coord for later use
    gap = alphabet('-')

    for i in range(alignment.get_alignment_length()):
        # convert to alphabet coordinates
        counts = Counter(
            alphabet(alignment[j, i].upper())
            for j in range(nrow)
            if j != refidx
        )
        pos, insert, label = next(labels)
        gapc = 0 if gap not in counts else counts[gap]
        maxc, minc, total = maxminsum(v for k, v in counts.items() if k != gap)
        # if everything is a gap (ie, total == 0), skip it
        if total == 0:
            yield None
        yield PosData(counts, gapc, insert, label, maxc, minc, pos, total)


def __poslabels(alignment, alphabet, refidx):
    refseq = str(alignment[refidx].seq)

    pos, insert = 0, 0
    for i, char in enumerate(refseq):
        if char not in Alphabet.GAP_CHARS:
            pos += 1
            insert = 0
        else:
            insert += 1
        yield (pos, insert, '%s%d' % (char.upper() if insert == 0 else '', pos))


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
