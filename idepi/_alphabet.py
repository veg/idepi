
from collections import defaultdict
from copy import deepcopy
from itertools import repeat

from numpy import dtype


__all__ = ['Alphabet']


class Alphabet(object):
    AMINO, DNA, STANFEL, CUSTOM = 0, 1, 2, 3

    # underscore must go first or re.compile blows up
    SPACE = '_.-='
    
    __DNA_ALPH   = 'ACGTUN-'
    __AMINO_ALPH = 'ACGILMPSTVDENQFWYHKRX-'
    __STANFEL_ALPH = {
        'A': 0, 'C': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0, 'P': 0, 'S': 0, 'T': 0, 'V': 0,
        'D': 1, 'E': 1, 'N': 1, 'Q': 1,
        'F': 2, 'W': 2, 'Y': 2,
        'H': 3, 'K': 3, 'R': 3,
        'X': 4,
        '-': 5 # this is a space, but we cant have 'keyword be an expression'
    }

    def __init__(self, mode=None, chars=None):
    
        def constant_factory(value):
            return repeat(value).next

        if mode is None:
           mode = self.AMINO

        if mode == Alphabet.STANFEL:
            d = Alphabet.__STANFEL_ALPH
            self.__list = []
            for i in xrange(max(d.values())):
                self.__list.append('[%s]' % ''.join(sorted([k if k not in Alphabet.SPACE else '' for k, v_ in d.items() if i == v_])))
        
        elif mode in (Alphabet.AMINO, Alphabet.DNA):
            alph = Alphabet.__AMINO_ALPH if mode == self.AMINO else Alphabet.__DNA_ALPH 
            d, self.__list = Alphabet.__dict_and_list(alph)

        elif mode == Alphabet.CUSTOM:
            if chars is None:
                raise ValueError('Custom Alphabet requires a `chars\' parameter')
            d, self.__list = Alphabet.__dict_and_list(chars)

        else:
            raise ValueError('mode must be one of Alphabet.AMINO, Alphabet.DNA, Alphabet.STANFEL, or Alphabet.CUSTOM')
        
        self.__dict = defaultdict(constant_factory(d['X']))
        self.__dict.update(d)

    def __len__(self):
        return len(self.__list)

    def __get__(self, idx):
        if type(idx) is int:
            return self.__list[idx]
        elif type(idx) is str:
            return self.__dict[idx]

    @staticmethod
    def __dict_and_list(alphabet):
        alph = alphabet.upper()
        d = dict([(alph[i], i) for i in xrange(len(alph))])
        l = [alph[i] if alph[i] not in Alphabet.SPACE else '[]' for i in xrange(len(alph))]
        return d, l

    def todict(self):
        return deepcopy(self.__dict)

    def tolist(self):
        return deepcopy(self.__list)

    def todtype(self):
        return dtype(('b1', (len(self.__list),)))
