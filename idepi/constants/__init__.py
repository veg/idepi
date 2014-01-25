
from __future__ import division, print_function

from Bio.Alphabet import Gapped
from Bio.Alphabet.IUPAC import extended_dna, extended_protein


GAPS = '_.-='

AminoAlphabet = Gapped(extended_protein)
DNAAlphabet = Gapped(extended_dna)
