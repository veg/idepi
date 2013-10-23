
from numpy import zeros

from sklearn.base import BaseEstimator, TransformerMixin

from idepi.filters import null_filter
from idepi.labeledmsa import LabeledMSA


__all__ = ['MSAVectorizer']


class MSAVectorizer(BaseEstimator, TransformerMixin):

    def __init__(self, encoder, filter=null_filter):
        self.__alignment_length = 0
        self.encoder = encoder
        self.filter = filter
        self.feature_names_ = []
        self.vocabulary_ = {}

    def fit(self, alignment):
        if not isinstance(alignment, LabeledMSA):
            raise ValueError("MSAVectorizers require a LabeledMSA")

        column_labels = list(alignment.labels)
        vocab = {}
        feature_names = []

        k = 0
        for i in range(alignment.get_alignment_length()):
            letters = self.filter(alignment[:, i])
            for ltr in letters:
                try:
                    j = self.encoder(ltr)
                    feature_name = '{0:s}{1:s}'.format(column_labels[i], self.encoder[j])
                    vocab[(i, j)] = k
                    feature_names.append(feature_name)
                    k += 1
                except KeyError:
                    pass

        self.__alignment_length = alignment.get_alignment_length()
        self.feature_names_ = feature_names
        self.vocabulary_ = vocab

        return self

    def transform(self, alignment):
        ncol = alignment.get_alignment_length()

        if ncol != self.__alignment_length:
            msg = 'alignment length ({0:d}) does not match the learned length ({1:d})'.format(
                ncol,
                self.__alignment_length)
            raise ValueError(msg)

        vocab = self.vocabulary_
        data = zeros((len(alignment), len(vocab)), dtype=int)

        if len(vocab) == 0:
            return data

        for i, seq in enumerate(alignment):
            seq_ = ''.join(ltr.upper() for ltr in str(seq.seq))
            for j, ltr in enumerate(seq_):
                try:
                    data[i, vocab[(j, self.encoder(ltr))]] = 1
                except KeyError:
                    pass

        return data

    def get_feature_names(self):
        return self.feature_names_
