
from numpy import zeros

from sklearn.base import BaseEstimator, TransformerMixin

from idepi.filters import null_filter
from idepi.labeledmsa import LabeledMSA


__all__ = ['MSAVectorizerPairwise']


class MSAVectorizerPairwise(BaseEstimator, TransformerMixin):

    def __init__(self, encoder, filter=null_filter, radius=0):
        if not isinstance(radius, int) or radius < 0:
            raise ValueError('radius expects a positive integer')
        self.__alignment_length = 0
        self.encoder = encoder
        self.filter = filter
        self.radius = radius
        self.feature_names_ = []
        self.vocabulary_ = {}

    def fit(self, alignment):
        if not isinstance(alignment, LabeledMSA):
            raise ValueError("MSAVectorizers require a LabeledMSA")

        valid_columns = [
            len(self.__filter(alignment[:, i])) > 0
            for i in range(alignment.get_alignment_length())]

        calls = set()

        for seq in alignment:
            seq_ = str(seq.seq).upper()
            for i, ltr1 in enumerate(seq_[:-1]):
                if not valid_columns[i]:
                    continue
                a = i + 1
                b = i + self.radius + 1
                try:
                    u = self.encoder(ltr1)
                    for j, ltr2 in enumerate(seq_[a:b], start=a):
                        if not valid_columns[j]:
                            continue
                        try:
                            v = self.encoder(ltr2)
                            calls.add((i, u, j, v))
                        except KeyError:
                            pass
                except KeyError:
                    pass

        column_labels = list(alignment.labels)
        feature_names = []
        vocab = {}

        for i, k in enumerate(sorted(calls)):
            vocab[k] = i
            idx1, ltr1, idx2, ltr2 = k
            feature_names.append('{0:s}{1:s}+{2:s}{3:s}'.format(
                column_labels[idx1],
                self.encoder[ltr1],
                column_labels[idx2],
                self.encoder[ltr2]
                ))

        self.__alignment_length = alignment.get_alignment_length()
        self.feature_names_ = feature_names
        self.vocabulary_ = vocab

        return self

    def transform(self, alignment):
        ncol = alignment.get_alignment_length()

        if ncol != self.__alignment_length:
            msg = 'alignment length ({0:d}) does not match the learned length ({1:d})'.format(
                ncol,
                self.__alignment_length
            )
            raise ValueError(msg)

        vocab = self.vocabulary_
        data = zeros((len(alignment), len(vocab)), dtype=int)

        if len(vocab) == 0:
            return data

        for row, seq in enumerate(alignment):
            seq_ = str(seq.seq).upper()
            for i, ltr1 in enumerate(seq_[:-1]):
                a = i + 1
                b = i + self.radius + 1
                try:
                    u = self.encoder(ltr1)
                    for j, ltr2 in enumerate(seq_[a:b], start=a):
                        try:
                            v = self.encoder(ltr2)
                            k = (i, u, j, v)
                            data[row, vocab[k]] = 1
                        except KeyError:
                            pass
                except KeyError:
                    pass

        return data

    def get_feature_names(self):
        return self.feature_names_
