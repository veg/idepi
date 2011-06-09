
from copy import deepcopy

from _seqtable import SeqTable


__all__ = ['Alphabet']


class Alphabet(object):
    AMINO, DNA, STANFEL = 0, 1, 2

    __STANFEL_ALPH = {
        'A': 0, 'C': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0, 'P': 0, 'S': 0, 'T': 0, 'V': 0,
        'D': 1, 'E': 1, 'N': 1, 'Q': 1,
        'F': 2, 'W': 2, 'Y': 2,
        'H': 3, 'K': 3, 'R': 3,
        'X': 4,
        '-': 5 # this is a space, but we cant have 'keyword be an expression'
    }

    def __init__(self, mode=None):
        if mode is None:
           mode = self.AMINO
        
        if mode not in (self.STANFEL, self.DNA, self.AMINO):
            raise ValueError('mode must be one of Alphabet.AMINO, Alphabet.DNA, or Alphabet.STANFEL')

        if mode == self.STANFEL:
            self.__dict = self.__STANFEL_ALPH
            self.__names = []
            for v in set(self.__dict.values()):
                self.__names.append('[%s]' % ''.join(sorted([k if k not in SeqTable.SPACE else '' for k, v_ in self.__dict.items() if v == v_])))
        else:
            alph = SeqTable.AMINO_ALPHABET if mode == self.AMINO else SeqTable.DNA_ALPHABET 
            self.__dict = dict([(alph[i], i) for i in xrange(len(alph))])
            self.__names = [alph[i] if alph[i] not in SeqTable.SPACE else '[]' for i in xrange(len(alph))]

    def todict(self):
        return deepcopy(self.__dict)

    def names(self):
        return deepcopy(self.__names)
