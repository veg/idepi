#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from sys import stderr
from Bio.Seq import Seq
from Bio.Alphabet import generic_nucleotide, generic_protein


__all__ = ['OrfList']


class OrfList(object):

    def __init__(self, seq):
        # accumulate all the potential start codon indices
        start_indices = []
        for codon in ("ATG", "AUG"):
            seq_ = seq.upper()
            prefix = 0
            while True:
                idx = seq_.find(codon)
                if idx >= 0:
                    start_indices.append(idx+prefix)
                    idx += 3
                    seq_ = seq_[idx:]
                    prefix += idx
                else:
                    break
        start_indices = sorted(start_indices)

        putative_seqs = []

        # find only those that terminate.. and truncate to its stop codon
        for idx in start_indices:
            p_dna_seq = seq[idx:].upper()
            p_dna_seq_ = Seq(p_dna_seq, generic_nucleotide)
            p_amino_seq = p_dna_seq_.translate()
            ter_idx = p_amino_seq.find('*')
            # if we don't find a termination codon, go ahead and add the full thing to the list
            if ter_idx < 0:
                putative_seqs.append((p_dna_seq_, p_amino_seq))
                continue
            # cut out small peptides (those <10 amino acids long)
            elif ter_idx < 30:
                continue
            # it remains, keep it
            p_seq = Seq(p_dna_seq[0:((ter_idx+1)*3)].upper(), generic_nucleotide)
            # truncate the stop "*" from the protein seq
            putative_seqs.append((p_seq, p_seq.translate()[:-1]))

        self.ORFs = sorted(putative_seqs, key=lambda d: len(d[1]), reverse=True)

    def __getitem__(self, key):
        return self.ORFs[key]

    def __len__(self):
        return len(self.ORFs)

    def __contains__(self, key):
        return True if key >= 0 and key < len(self.ORFs) else False
