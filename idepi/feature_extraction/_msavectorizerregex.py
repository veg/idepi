
from numpy import zeros

from sklearn.base import BaseEstimator, TransformerMixin

from idepi.constants import GAPS
from idepi.labeledmsa import LabeledMSA


__all__ = ['MSAVectorizerRegex']


class MSAVectorizerRegex(BaseEstimator, TransformerMixin):

    def __init__(self, regex, regex_length=-1, name=''):
        self.name = name
        self.regex = regex
        self.regex_length = regex_length
        self.feature_names_ = []
        self.vocabulary_ = {}

    def fit(self, alignment):
        if not isinstance(alignment, LabeledMSA):
            raise ValueError("MSAVectorizers require a LabeledMSA")

        calls = set()

        for idx, seq in enumerate(alignment):
            # do NOT convert to encoder coords
            cols = []
            ltrs = []
            for col, ltr in enumerate(str(seq.seq)):
                if ltr not in GAPS:
                    cols.append(col)
                    ltrs.append(ltr)

            seq_ = ''.join(ltrs)

            for m in self.regex.finditer(seq_):
                start = m.start(0)
                idx = cols[start]
                calls.add(idx)
                if self.regex_length >= 0 and len(m.group(0)) != self.regex_length:
                    raise ValueError(
                        "supplied regex_length incorrect for: '{0}'".format(m.group(0))
                        )

        column_labels = list(alignment.labels)
        feature_names = []
        vocab = {}

        for j, idx in enumerate(sorted(calls)):
            vocab[idx] = j
            feature_names.append('{0:s}({1:s})'.format(self.name, column_labels[idx]))

        self.feature_names_ = feature_names
        self.vocabulary_ = vocab

        return self

    def transform(self, alignment):
        vocab = self.vocabulary_
        data = zeros((len(alignment), len(vocab)), dtype=int)

        if len(vocab) == 0:
            return data

        for i, seq in enumerate(alignment):
            # do NOT convert to encoder coords
            cols = []
            ltrs = []
            for col, ltr in enumerate(str(seq.seq)):
                if ltr not in GAPS:
                    cols.append(col)
                    ltrs.append(ltr)

            seq_ = ''.join(ltrs)

            # generate a list of all aln col idx at which the pattern is found,
            # the sorted should be unnecessary (finditer scans l-to-r), but I'm paranoid
            matches = sorted(cols[m.start(0)] for m in self.regex.finditer(seq_))

            # match all idx pairs to the ones in filtercalls,
            # and set the data matrix appropriately
            for idx in matches:
                try:
                    data[i, vocab[idx]] = 1
                except KeyError:
                    pass

        return data

    def get_feature_names(self):
        return self.feature_names_
