
from numpy import zeros

from ._filters import nofilter
from ._posstream import posstream
from .._common import base_10_to_n, base_26_to_alph


__all__ = ['DataBuilder1D']


class DataBuilder1D(object):

    def __init__(self, alignment, alphabet, refidx, filter=nofilter):
        self.__alphabet = alphabet
        self.__filtercalls = []
        self.__labels = []

        # evaluate each position in the stream and generate the column labels
        for p in posstream(alignment, alphabet, refidx):
            chars = filter(p)
            self.__filtercalls.append(chars)
            for char in chars:
                insert = base_26_to_alph(base_10_to_n(p.insert, 26))
                self.__labels.append('%s%s%s' % (p.label, insert, char))

        self.__length = sum(len(chars) for chars in self.__filtercalls)

        assert self.__length == len(self.__labels)

    def __len__(self):
        return self.__length

    def __call__(self, alignment, refidx=None):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        if alignment.get_alignment_length() != len(self.__filtercalls):
            msg = 'alignment length (%d) does not match the learned length (%d)' % (
                alignment.get_alignment_length(),
                len(self.__filtercalls)
            )
            raise ValueError(msg)

        data = zeros((
            len(alignment) - (0 if refidx is None else 1),
            len(self)
        ), dtype=int)

        col = 0
        for i, chars in enumerate(self.__filtercalls):
            for j, char in enumerate(chars, start=col):
                # handle the refidx case
                column = (alignment[k, i] for k in range(len(alignment)) if k != refidx)
                for k, c in enumerate(column):
                    data[k, j] = c == char
            col += len(chars)

        return data

    @property
    def labels(self):
        return self.__labels
