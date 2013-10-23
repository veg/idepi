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

from __future__ import division, print_function

from os import environ
from os.path import exists, join
from re import compile as re_compile
from subprocess import Popen, PIPE

from Bio.Seq import Seq


__all__ = ['HMMER']


class HMMER:
    AMINO, DNA, RNA = 0, 1, 2
    FASTA, EMBL, GENBANK, UNIPROT = 'FASTA', 'EMBL', 'Genbank', 'Uniprot'
    STOCKHOLM, PFAM, A2M, PSIBLAST = 'Stockholm', 'Pfam', 'A2M', 'PSIBLAST'

    def __init__(self, alignbin='hmmalign', buildbin='hmmbuild'):
        self.__alignbin = alignbin
        self.__buildbin = buildbin

        for bin in (self.__alignbin, self.__buildbin):
            ex = False
            if exists(self.__alignbin):
                ex = True
            else:
                for path in environ['PATH'].split(':'):
                    if exists(join(path, bin)):
                        ex = True
            if not ex:
                raise Exception('Executable %s not found.' % bin)

    def align(self, hmmfile, seqfile, output=None, mapali=None, trim=False, alphabet=None, informat=None, outformat=None):

        if alphabet is None:
            alphabet = HMMER.AMINO
        if informat is None:
            informat = HMMER.FASTA
        if outformat is None:
            outformat = HMMER.STOCKHOLM

        if alphabet not in (HMMER.AMINO, HMMER.DNA, HMMER.RNA):
            raise ValueError('alphabet needs to be one of the idepi-provided constants: AMINO, DNA, RNA')
        if informat not in (HMMER.FASTA, HMMER.EMBL, HMMER.GENBANK, HMMER.UNIPROT):
            raise ValueError('informat needs to be one of the idepi-provided constants: FASTA, EMBL, GENBANK, UNIPROT')
        if outformat not in (HMMER.STOCKHOLM, HMMER.PFAM, HMMER.A2M, HMMER.PSIBLAST):
            raise ValueError('outformat needs to be one of the idepi-provided constants: STOCKHOLM, PFAM, A2M, PSIBLAST')

        # I don't think this is necessary code
#         tmp = False
#         if output is None:
#             fd, output = mkstemp(); close(fd)
#             tmp = True # TODO something with this variable...?

        args = [self.__alignbin]
        if mapali is not None:
            args.extend(['--mapali', mapali])
        if trim is True:
            args.append('--trim')
        if alphabet is not None:
            if alphabet == HMMER.AMINO:
                args.append('--amino')
            elif alphabet == HMMER.DNA:
                args.append('--dna')
            elif alphabet == HMMER.RNA:
                args.append('--rna')
        if informat is not None:
            args.extend(['--informat', informat])
        if outformat is not None:
            args.extend(['--outformat', outformat])
        if output is not None:
            args.extend(['-o', output])
        args.extend([hmmfile, seqfile])

        kwargs = {}
        if output is None:
            kwargs['stdout'] = PIPE

        proc = Popen(args, close_fds=True, stderr=PIPE, **kwargs)
        out, err = proc.communicate()

        if proc.returncode is not 0:
            return RuntimeError(err)

        if output is None:
            return out

    def build(self, hmmfile, alignmentfile, name=None, logfile=None, annotated=None, alphabet=None):

        args = [self.__buildbin]
        if name is not None:
            args.extend(['-n', name])
        if logfile is not None:
            args.extend(['-o', logfile])
        if annotated is not None:
            args.extend(['-O', annotated])
        if alphabet is not None:
            if alphabet == HMMER.AMINO:
                args.append('--amino')
            elif alphabet == HMMER.DNA:
                args.append('--dna')
            elif alphabet == HMMER.RNA:
                args.append('--rna')
            else:
                raise ValueError('unknown alphabet passed to HMMER')
        args.extend([hmmfile, alignmentfile])

        kwargs = {}
        if logfile is None:
            kwargs['stdout'] = PIPE

        proc = Popen(args, close_fds=True, stderr=PIPE, **kwargs)
        out, err = proc.communicate()

        if proc.returncode is not 0:
            return RuntimeError(err)

        if logfile is None:
            return out

    @staticmethod
    def valid(record, is_dna=False):
        if is_dna:
            regexp = re_compile(r'[^ACGT]')
        else:
            regexp = re_compile(r'[^ACDEFGHIKLMNPQRSTVWY]')

        seq = regexp.sub('', str(record.seq))

        record.letter_annotations.clear()
        record.seq = Seq(seq, record.seq.alphabet)

        return record
