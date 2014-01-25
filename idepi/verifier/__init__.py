
from __future__ import division, print_function

from re import compile as re_compile, I as re_I


__all__ = [
    'VerifyError',
    'Verifier'
    ]


def verify_alphabet(sequence, alphabet=None):
    if alphabet is None:
        alphabet = sequence.alphabet
    letters = alphabet.letters
    if not letters:
        return True
    for letter in str(sequence):
        if letter not in letters:
            return False
    return True


class VerifyError(Exception):
    pass


class Verifier:

    def __init__(self, source, alphabet):
        self.set_alphabet(alphabet)
        self.__records = []
        self.__source = enumerate(source)

    def __iter__(self):
        return self()

    def __catchup(self):
        for i, record in self.__records:
            record.seq.alphabet = self.__alphabet
            if self.__regexp is not None and self.__regexp.match(str(record.seq)):
                raise VerifyError("invalid alphabet for sequence {0:d}".format(i))
        while self.__records:
            _, record = self.__records.pop(0)
            yield record

    def __call__(self):
        for item in self.__source:
            self.__records.append(item)
            for record in self.__catchup():
                yield record
        for record in self.__catchup():
            yield record

    def set_alphabet(self, alphabet):
        self.__alphabet = alphabet
        if alphabet.letters:
            self.__regexp = re_compile(r'[^{0:s}]'.format(alphabet.letters), re_I)
        else:
            self.__regexp = None
