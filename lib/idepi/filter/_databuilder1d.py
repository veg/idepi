
from numpy import zeros

from ._posstream import posstream
from .._common import base_10_to_n, base_26_to_alph

class DataBuilder1D(object):

    def __init__(self):
        self.__alphabet = None
        self.__filtercalls = []
        self.__labels = []
        self.__length = None

    def __learn(self, alignment, alphabet, filter, refidx, radius):

        self.__alphabet = alphabet

        # evaluate each position in the stream and generate the column labels
        for p in posstream(alignment, alphabet, refidx):
            chars = filter(p)
            self.__filtercalls.append(chars)
            for char in chars:
                insert = base_26_to_alph(base_10_to_n(p.insert, 26))
                self.__labels.append('%s%s%s' % (p.label, insert, char))

        self.__length = sum(len(chars) for chars in self.__filtercalls)

    def __len__(self):
        return self.__length

    def filter(self, alignment, refidx=None):
        if self.__length is None:
            raise RuntimeError('no filter model computed! call learn()')

        if alignment.get_alignment_length() != len(self.__filtercalls):
            msg = 'alignment length (%d) does not match the learned length (%d)' % (
                alignment.get_alignment_length(),
                len(self.__filtercalls)
            )
            raise ValueError(msg)

        data = zeros((len(alignment), len(self)), dtype=int)

        col = 0
        for i, chars in enumerate(self.__filtercalls):
            for j, char in enumerate(chars, start=col):
                for k, c in enumerate(alignment[:, i]):
                    data[k, j] = c == char
            col += len(chars)

        return data

    @property
    def labels(self):
        return self.__labels

    def learn(self, alignment, alphabet, filter, refidx, radius=0):
        self.__compute(alignment, alphabet, filter, refidx, radius)
        return self.filter(alignment)
