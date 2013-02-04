
from numpy import zeros

from BioExt.collections import OrderedDict

from ._posstream import posstream
from ..alphabet import Alphabet


__all__ = ['DataBuilderRegex']


class DataBuilderRegex:

    def __init__(self, alignment, alphabet, refidx, regex, regex_length=-1, label=''):

        self.__alphabet = alphabet
        self.__labels = []
        self.__regex = regex
        self.__regex_length = regex_length

        # save these around
        pstream = list(posstream(alignment, alphabet, refidx))

        calls = set()

        for i, seq in enumerate(alignment):
            if i == refidx:
                continue

            # do NOT convert to alphabet coords
            cols, chars = zip(*[
                (col, char)
                for col, char in enumerate(seq)
                if char not in Alphabet.GAPS
            ])

            seqp = ''.join(chars)

            for m in self.__regex.finditer(seqp):
                start = m.start(0)
                idx = cols[start]
                calls.add(idx)
                if regex_length >= 0 and len(m.group(0)) != regex_length:
                    raise ValueError(
                        "supplied regex_length incorrect for: '{0}'".format(m.group(0))
                        )

        self.__filtercalls = OrderedDict(
            (v, k)
            for k, v in enumerate(sorted(calls))
            )

        for idx in sorted(self.__filtercalls.keys()):
            self.__labels.append('%s(%s)' % (
                label, pstream[idx].label
            ))

        self.__length = len(self.__filtercalls)

        assert self.__length == len(self.__labels)

    def __len__(self):
        return self.__length

    def __call__(self, alignment, refidx=None, globber=None, normalize=False):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        if normalize and self.__regex_length < 0:
            raise ValueError(
                'normalize requires regex_length to be provided during initialization'
                )

        if globber is None:
            nrow = len(alignment) - (0 if refidx is None else 1)
        else:
            nrow = len(globber)

        data = zeros(
            (nrow, len(self)),
            dtype=int if globber is None else float
            )

        if len(self) == 0:
            return data

        alignment_ = (seq for i, seq in enumerate(alignment) if i != refidx)

        if normalize:
            coverage = zeros((len(self),), dtype=int)
            gaps = set(Alphabet.GAPS)

        for i, seq in enumerate(alignment_):
            if globber is None:
                r, weight = i, 1
            else:
                r, weight = globber[seq.id]

            # do NOT convert to alphabet coords
            cols, chars = zip(*[
                (col, char)
                for col, char in enumerate(seq)
                if char not in Alphabet.GAPS
            ])

            seq_ = ''.join(chars)

            if normalize:
                for idx, j in self.__filtercalls.items():
                    lwr, upr = idx, idx + self.__regex_length
                    if not set(seq[lwr:upr]) & gaps:
                        coverage[j] += 1

            for m in self.__regex.finditer(seq_):
                start = m.start(0)
                idx = cols[start]
                if idx in self.__filtercalls:
                    j = self.__filtercalls[idx]
                    data[r, j] += weight

        if normalize:
            return data / coverage

        return data

    @property
    def labels(self):
        return self.__labels
