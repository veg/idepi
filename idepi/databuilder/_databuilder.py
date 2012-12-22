
from numpy import zeros

from ._posstream import posstream
from ..filter import nofilter


__all__ = ['DataBuilder']


class DataBuilder:

    def __init__(self, alignment, alphabet, refidx, filter=nofilter):
        self.__alphabet = alphabet
        self.__filtercalls = []
        self.__labels = []

        # evaluate each position in the stream and generate the column labels,
        # being careful to use the alphabet repr
        for p in posstream(alignment, alphabet, refidx):
            chars = filter(p)
            self.__filtercalls.append(chars)
            for char in chars:
                self.__labels.append('%s%s' % (p.label, alphabet[char]))

        self.__length = sum(len(chars) for chars in self.__filtercalls)

        assert self.__length == len(self.__labels)

    def __len__(self):
        return self.__length

    def __call__(self, alignment, refidx=None, globber=None):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        if alignment.get_alignment_length() != len(self.__filtercalls):
            msg = 'alignment length (%d) does not match the learned length (%d)' % (
                alignment.get_alignment_length(),
                len(self.__filtercalls)
            )
            raise ValueError(msg)

        if globber is None:
            nrow = len(alignment) - (0 if refidx is None else 1)
        else:
            nrow = len(globber)

        data = zeros((nrow, len(self)), dtype=int)

        col = 0
        for i, chars in enumerate(self.__filtercalls):
            for j, char in enumerate(chars, start=col):
                # handle the refidx case,
                # and convert to alphabet coordinates
                if globber is None:
                    column = enumerate(
                        self.__alphabet(alignment[k, i])
                        for k in range(len(alignment))
                        if k != refidx
                        )
                else:
                    column = (
                        (globber[alignment[k].id], self.__alphabet(alignment[k, i]))
                        for k in range(len(alignment))
                        if k != refidx
                        )
                for k, c in column:
                    data[k, j] += c == char
            col += len(chars)

        return data

    @property
    def labels(self):
        return self.__labels
