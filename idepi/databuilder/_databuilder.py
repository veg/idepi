
from numpy import zeros

from ..alphabet import Alphabet
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

    def __call__(self, alignment, refidx=None, globber=None, normalize=False):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        ncol = alignment.get_alignment_length()

        if ncol != len(self.__filtercalls):
            msg = 'alignment length (%d) does not match the learned length (%d)' % (
                ncol,
                len(self.__filtercalls)
            )
            raise ValueError(msg)

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
            for j, chars in enumerate(self.__filtercalls):
                if normalize and seq_[j] not in Alphabet.GAPS:
                    lwr, upr = col, col + len(chars)
                    coverage[lwr:upr] += 1
                for k, char in enumerate(chars, start=col):
                    # convert to alphabet coordinates
                    data[r, k] += self.__alphabet(seq_[j]) == char
                col += len(chars)

        if normalize:
            return data / coverage

        return data

    @property
    def labels(self):
        return self.__labels
