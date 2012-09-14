
from numpy import zeros

from ._posstream import posstream
from ..alphabet import Alphabet


__all__ = ['DataBuilderRegex']


class DataBuilderRegex(object):

    def __init__(self, alignment, alphabet, refidx, regex, label=''):

        self.__alphabet = alphabet
        self.__labels = []
        self.__regex = regex

        # save these around
        pstream = list(posstream(alignment, alphabet, refidx))

        calls = set()

        for idx, seq in enumerate(alignment):
            if idx == refidx:
                continue

            # do NOT convert to alphabet coords
            cols, chars = zip(*[
                (col, char)
                for col, char in enumerate(seq)
                if char not in Alphabet.GAP_CHARS
            ])

            seqp = ''.join(chars)

            for m in self.__regex.finditer(seqp):
                start = m.start(0)
                idx = cols[start]
                calls.add(idx)

        self.__filtercalls = dict((v, k) for k, v in enumerate(sorted(calls)))

        for idx in sorted(self.__filtercalls.keys()):
            self.__labels.append('%s(%s)' % (
                label, pstream[idx].label
            ))

        self.__length = len(self.__filtercalls)

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

        alignmentp = (seq for i, seq in enumerate(alignment) if i != refidx)

        for i, seq in enumerate(alignmentp):
            # do NOT convert to alphabet coords
            cols, chars = zip(*[
                (col, char)
                for col, char in enumerate(seq)
                if char not in Alphabet.GAP_CHARS
            ])

            seqp = ''.join(chars)

            for m in self.__regex.finditer(seqp):
                start = m.start(0)
                idx = cols[start]
                if idx in self.__filtercalls:
                    j = self.__filtercalls[idx]
                    data[i, j] = True

        return data

    @property
    def labels(self):
        return self.__labels
