
from numpy import zeros

from ._posstream import posstream
from ..alphabet import Alphabet
from ..filter import nofilter


__all__ = ['DataBuilderPairwise']


def clamp(min, max, val):
    if val >= max:
        return max - 1
    elif val < min:
        return min
    else:
        return val


class DataBuilderPairwise:

    def __init__(self, alignment, alphabet, refidx, filter=nofilter, radius=0):
        if radius < 0:
            raise ValueError('radius expects a positive integer')

        self.__alphabet = alphabet
        self.__labels = []

        if radius == 0:
            self.__filtercalls = []
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
                if char not in Alphabet.GAPS
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
                    elif pidx1 > pidx2:
                        key = (pidx2, pidx1)
                        val = (char2, char1)
                    else:
                        continue
                    if key not in calls:
                        calls[key] = set()
                    calls[key].add(val)

        # this correctly sorts the character pairs
        def alphkey(x):
            return len(alphabet) * x[0] + x[1]

        # and this correctly sorts the idx pairs
        def idxkey(x):
            return len(pstream) * x[0][0] + x[0][1]

        # turn the dict of sets into a sorted list of sorted lists
        self.__filtercalls = sorted(
            (
                (k, sorted(v, key=alphkey))
                for k, v in calls.items()
                ),
            key=idxkey
            )

        for ab, pairs in self.__filtercalls:
            a, b = ab
            for pair in pairs:
                self.__labels.append(
                    '%s%s+%s%s' % (
                        # convert to alphabet repr
                        pstream[a].label, alphabet[pair[0]],
                        pstream[b].label, alphabet[pair[1]]
                        )
                    )

        self.__length = sum(len(pairs) for _, pairs in self.__filtercalls)

        assert self.__length == len(self.__labels)

    def __len__(self):
        return self.__length

    def __call__(self, alignment, refidx=None, globber=None, normalize=False):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        if globber is None:
            nrow = len(alignment) - (0 if refidx is None else 1)
        else:
            nrow = len(globber)

        data = zeros((nrow, len(self)), dtype=int)

        if len(self) == 0:
            return data

        alignment_ = (seq for i, seq in enumerate(alignment) if i != refidx)

        if normalize:
            coverage = zeros((len(self),), dtype=int)

        for i, seq in enumerate(alignment_):
            if globber is None:
                r = i
            else:
                r = globber[seq.id]
            seq_ = ''.join(c.upper() for c in seq)
            col = 0
            for ab, pairs in self.__filtercalls:
                a, b = ab
                if (normalize and
                    seq_[a] not in Alphabet.GAPS and
                    seq_[b] not in Alphabet.GAPS):
                    lwr, upr = col, col + len(pairs)
                    coverage[lwr:upr] += 1
                for j, pair in enumerate(pairs, start=col):
                    c, d = pair
                    if (self.__alphabet(seq_[a]) == c and
                        self.__alphabet(seq_[b]) == d):
                        data[r, j] += True
                col += len(pairs)

        if normalize:
            return data / coverage

        return data

    @property
    def labels(self):
        return self.__labels
