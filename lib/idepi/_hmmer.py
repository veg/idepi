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

import logging

from os import close, environ
from os.path import exists, join
from subprocess import Popen, PIPE
from tempfile import mkstemp


__all__ = ['Hmmer']


class Hmmer(object):
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

    def align(self, hmmfile, seqfile, output=None, allcol=False, mapali=None, trim=False, alphabet=None, informat=None, outformat=None):

        if alphabet is None:
            alphabet = Hmmer.AMINO
        if informat is None:
            informat = Hmmer.FASTA
        if outformat is None:
            outformat = Hmmer.STOCKHOLM

        if alphabet not in (Hmmer.AMINO, Hmmer.DNA, Hmmer.RNA):
            raise ValueError('alphabet needs to be one of the idepi-provided constants: AMINO, DNA, RNA')
        if informat not in (Hmmer.FASTA, Hmmer.EMBL, Hmmer.GENBANK, Hmmer.UNIPROT):
            raise ValueError('informat needs to be one of the idepi-provided constants: FASTA, EMBL, GENBANK, UNIPROT')
        if outformat not in (Hmmer.STOCKHOLM, Hmmer.PFAM, Hmmer.A2M, Hmmer.PSIBLAST):
            raise ValueError('outformat needs to be one of the idepi-provided constants: STOCKHOLM, PFAM, A2M, PSIBLAST')

        tmp = False
        if output is None:
            fd, output = mkstemp(); close(fd)
            tmp = True

        args = [self.__alignbin]
        if allcol is True:
            args.append('--allcol')
        if mapali is not None:
            args.extend(['--mapali', mapali])
        if trim is True:
            args.append('--trim')
        if alphabet is not None:
            if alphabet == Hmmer.AMINO:
                args.append('--amino')
            if alphabet == Hmmer.DNA:
                args.append('--dna')
            if alphabet == Hmmer.RNA:
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

    def build(self, hmmfile, alignmentfile, name=None, logfile=None, annotated=None):

        args = [self.__buildbin]
        if name is not None:
            args.extend(['-n', name])
        if logfile is not None:
            args.extend(['-o', logfile])
        if annotated is not None:
            args.extend(['-O', annotated])
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
