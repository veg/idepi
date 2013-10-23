
from Bio.Align import MultipleSeqAlignment

from BioExt.collections import Counter

from idepi.constants import GAPS
from idepi._common import base_10_to_n, base_26_to_alph


__all__ = [
    'LabeledMSA',
    'column_labels'
    ]


def column_labels(alignment, refidx):
    refseq = str(alignment[refidx].seq)
    pos, insert = 0, 0
    for i, char in enumerate(refseq):
        if char not in GAPS:
            pos += 1
            insert = 0
        else:
            insert += 1
        ins = base_26_to_alph(base_10_to_n(insert, 26))
        yield pos, '{0:s}{1:d}{2:s}'.format(char.upper() if insert == 0 else '', pos, ins)


class LabeledMSA(MultipleSeqAlignment):

    @staticmethod
    def from_msa_with_ref(msa, refidx):
        positions, labels = [], []
        for pos, label in column_labels(msa, refidx):
            positions.append(pos)
            labels.append(label)
        msa2 = msa[:refidx]
        msa2.extend(msa[refidx + 1:])
        return LabeledMSA(
            msa2,
            labels,
            positions
            )

    def __init__(self, msa, labels, positions):
        if not isinstance(msa, MultipleSeqAlignment):
            raise TypeError("invalid msa type")
        ncol = msa.get_alignment_length()
        if ncol > 0 and (len(labels) != ncol or len(positions) != ncol):
            raise ValueError("all arguments must share the same column space")
        self.__labels = labels
        self.__positions = positions
        super(LabeledMSA, self).__init__(iter(msa), msa._alphabet)

    def __getitem__(self, index):
        if isinstance(index, (int, slice)):
            return LabeledMSA(
                super(LabeledMSA, self).__getitem__(index),
                self.__labels,
                self.__positions
                )
        elif len(index) != 2 or not all(isinstance(idx, (int, slice)) for idx in index):
            raise TypeError("invalid index type")

        _, col_index = index
        if isinstance(col_index, int):
            return super(LabeledMSA, self).__getitem__(index)
        elif isinstance(col_index, slice):
            return LabeledMSA(
                super(LabeledMSA, self).__getitem__(index),
                self.__labels[col_index],
                self.__positions[col_index]
                )
        else:
            raise TypeError("invalid index type")

    def __add__(self, other):
        if not isinstance(other, LabeledMSA):
            raise NotImplementedError
        return LabeledMSA(
            super(LabeledMSA, self).__add__(other),
            self.__labels + other.__labels,
            self.__positions + other.__positions
            )

    def get_alignment_length(self):
        return len(self.__labels)

    @property
    def labels(self):
        return iter(self.__labels)

    @property
    def positions(self):
        return iter(self.__positions)
