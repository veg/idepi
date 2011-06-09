
from Bio.SeqRecord import SeqRecord

from _orflist import OrfList


__all__ = ['AbRecord']


class AbRecord(object):

    def __init__(self, id, seq, subtype, ab, ic50):
        self.id = id
        self.dna_seq, self.amino_seq = OrfList(seq)[0]
        self.subtype = subtype.upper()
        self.antibody = ab
        ic50 = ic50.strip()
        if '>' in ic50:
            ic50 = 25.
        else:
            ic50 = float(ic50)
        self.ic50 = ic50

    def to_SeqRecord(self, dna=False):
        if dna:
            _seq = self.dna_seq
        else:
            _seq = self.amino_seq
        return SeqRecord(_seq, '%s|%s|%s|%s' % (self.id, self.subtype, self.antibody, self.ic50))
