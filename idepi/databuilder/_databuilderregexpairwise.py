
from operator import itemgetter

from numpy import zeros

from BioExt.collections import OrderedDict

from ._posstream import posstream
from ..alphabet import Alphabet


__all__ = ['DataBuilderRegexPairwise']


class DataBuilderRegexPairwise:

    def __init__(self, alignment, alphabet, regex, regex_length=-1, label=''):

        self.__alphabet = alphabet
        self.__labels = []
        self.__regex = regex
        self.__regex_length = regex_length

        # save these around
        pstream = list(posstream(alignment, alphabet))

        calls = set()

        for idx, seq in enumerate(alignment):
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

        self.__filtercalls = OrderedDict()

        # j is the column idx in the resultant data matrix
        j = 0
        for idx1 in sorted(calls):
            for idx2 in sorted(calls):
                if idx2 <= idx1:
                    continue
                self.__filtercalls[(idx1, idx2)] = j
                j += 1

        # sort on the value, which is the column idx (see above)
        for k, _ in sorted(self.__filtercalls.items(), key=itemgetter(1)):
            idx1, idx2 = k
            self.__labels.append('%s(%s+%s)' % (
                label, pstream[idx1].label, pstream[idx2].label
            ))

        self.__length = len(self.__filtercalls)

        assert self.__length == len(self.__labels)

    def __len__(self):
        return self.__length

    def __call__(self, alignment, globber=None, normalize=False):
        if self.__length is None:
            raise RuntimeError('no filter model computed! programmer error!')

        if normalize and self.__regex_length < 0:
            raise ValueError(
                'normalize requires regex_length to be provided during initialization'
                )

        if globber is None:
            nrow = len(alignment)
        else:
            nrow = len(globber)

        data = zeros(
            (nrow, len(self)),
            dtype=int if globber is None else float
            )

        if len(self) == 0:
            return data

        if normalize:
            coverage = zeros((len(self),), dtype=int)
            gaps = set(Alphabet.GAPS)

        for i, seq in enumerate(alignment):
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

            # generate a list of all aln col idx at which the pattern is found,
            # the sorted should be unnecessary (finditer scans l-to-r), but I'm paranoid
            matches = sorted(cols[m.start(0)] for m in self.__regex.finditer(seq_))

            if normalize:
                for item in self.__filtercalls.items():
                    idxs, j = item
                    lwr1, lwr2 = idxs
                    upr1 = lwr1 + self.__regex_length
                    upr2 = lwr2 + self.__regex_length
                    if not set(seq[lwr1:upr1]) & set(seq[lwr2:upr2]) & gaps:
                        coverage[j] += 1

            # match all idx pairs to the ones in filtercalls,
            # and set the data matrix appropriately
            for idx1 in matches:
                for idx2 in matches:
                    if idx2 <= idx1:
                        continue
                    k = (idx1, idx2)
                    if k in self.__filtercalls:
                        j = self.__filtercalls[k]
                        data[r, j] += weight

        if normalize:
            return data / coverage

        return data

    @property
    def labels(self):
        return self.__labels
