
from numpy import zeros

from ._posstream import posstream
from ..alphabet import Alphabet
from ..filter import nofilter
from .._common import base_10_to_n, base_26_to_alph


__all__ = ['DataBuilder2D']


def clamp(min, max, val):
    if val >= max:
        return max - 1
    elif val < min:
        return min
    else:
        return val


class DataBuilder2D(object):

    def __init__(self, alignment, alphabet, refidx, filter=nofilter, radius=0):
        if radius < 0:
            raise ValueError('radius expects a positive integer')

        self.__alphabet = alphabet
        self.__filtercalls = {}
        self.__labels = []

        if radius == 0:
            self.__length = 0
            return

        # save these around
        pstream = list(posstream(alignment, alphabet, refidx))

        calls = {}

        for idx, seq in enumerate(alignment):
            if idx == refidx:
                continue
            # take care to convert to alphabet coords
            seqp = [
                (col, alphabet(char))
                for col, char in enumerate(seq)
                if char not in Alphabet.GAP_CHARS
            ]
            for i, colchar in enumerate(seqp):
                pidx1, char1 = colchar
                # if the filter rejects this position, skip
                if not filter(pstream[pidx1]):
                    continue
                lwr = clamp(0, len(seqp), i - radius)
                upr = clamp(0, len(seqp), i + radius + 1)
                for pidx2, char2 in (seqp[k] for k in range(lwr, upr)):
                    # if the filter rejects this position, skip
                    if not filter(pstream[pidx2]):
                        continue
                    # the array should be upper-triangular
                    if pidx1 < pidx2:
                        key = (pidx1, pidx2)
                        val = (char1, char2)
                    elif pidx2 < pidx1:
                        key = (pidx2, pidx1)
                        val = (char2, char1)
                    else:
                        continue
                    if key not in calls:
                        calls[key] = set()
                    calls[key].add(val)

        # this correctly sorts the alphabet characters
        def keyfn1(x):
            return len(alphabet) * x[0] + x[1]

        # this performs a 2D sorting of the
        # list of position pairs
        def keyfn2(x):
            return len(pstream) * x[0][0] + x[0][1]

        # turn the dict of sets into a sorted list of sorted lists
        self.__filtercalls = sorted(((k, sorted(v, key=keyfn1)) for k, v in calls.items()), key=keyfn2)

        for ab, pairs in self.__filtercalls:
            a, b = ab
            for pair in pairs:
                ins1 = base_26_to_alph(base_10_to_n(pstream[a].insert, 26))
                ins2 = base_26_to_alph(base_10_to_n(pstream[b].insert, 26))
                self.__labels.append('%s%s%s+%s%s%s' % (
                    # convert to alphabet repr
                    pstream[a].label, ins1, alphabet[pair[0]],
                    pstream[b].label, ins2, alphabet[pair[1]]
                ))

        self.__length = sum(len(pairs) for _, pairs in self.__filtercalls)

        assert self.__length == len(self.__labels)

    def __len__(self):
        return self.__length

    def __call__(self, alignment, refidx=None):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        data = zeros((
            len(alignment) - (0 if refidx is None else 1),
            len(self)
        ), dtype=int)

        col = 0
        for ab, pairs in self.__filtercalls:
            a, b = ab
            for j, pair in enumerate(pairs, start=col):
                # handle the refidx case,
                # and convert to alphabet coordinates
                column = (
                    (
                        self.__alphabet(alignment[k, a]),
                        self.__alphabet(alignment[k, b])
                    )
                    for k in range(len(alignment))
                    if k != refidx
                )
                for k, c in enumerate(column):
                    data[k, j] = (c == pair)
            col += len(pairs)

        return data

    @property
    def labels(self):
        return self.__labels
