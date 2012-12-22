
from operator import itemgetter

from numpy import zeros

from ._posstream import posstream
from ..alphabet import Alphabet


__all__ = ['DataBuilderRegexTriplewise']


class DataBuilderRegexTriplewise:

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

        self.__filtercalls = dict()

        # j is the column idx in the resultant data matrix
        j = 0
        for idx1 in sorted(calls):
            for idx2 in sorted(calls):
                if idx2 <= idx1:
                    continue
                for idx3 in sorted(calls):
                    if idx3 <= idx1 or idx3 <= idx2:
                        continue
                    self.__filtercalls[(idx1, idx2, idx3)] = j
                    j += 1

        # sort on the value, which is the column idx (see above)
        for k, _ in sorted(self.__filtercalls.items(), key=itemgetter(1)):
            idx1, idx2 = k
            self.__labels.append('%s(%s+%s+%s)' % (
                label, pstream[idx1].label, pstream[idx2].label, pstream[idx3].label
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

            # generate a list of all aln col idx at which the pattern is found,
            # the sorted should be unnecessary (finditer scans l-to-r), but I'm paranoid
            matches = sorted(cols[m.start(0)] for m in self.__regex.finditer(seqp))

            # match all idx pairs to the ones in filtercalls,
            # and set the data matrix appropriately
            for idx1 in matches:
                for idx2 in matches:
                    if idx2 <= idx1:
                        continue
                    for idx3 in matches:
                        if idx3 <= idx1 or idx3 <= idx2:
                            continue
                        k = (idx1, idx2, idx3)
                        if k in self.__filtercalls:
                            j = self.__filtercalls[k]
                            data[i, j] = True

        return data

    @property
    def labels(self):
        return self.__labels
