#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# sto2fa.py :: converts a stockholm multiple sequence alignment file to fasta
# format.
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

from argparse import ArgumentParser, FileType
from sys import argv, exit, stdin, stdout

from Bio import AlignIO


def main(args=None):

    if args is None:
        args = argv[1:]

    parser = ArgumentParser(description='convert stockholm file to FASTA')
    parser.add_argument('STOCKHOLMFILE', type=FileType('r'))
    ns = parser.parse_args(args)

    alignments = AlignIO.parse(ns.STOCKHOLMFILE, 'stockholm')
    AlignIO.write(alignments, stdout, 'fasta')

    if ns.STOCKHOLMFILE != stdin:
        ns.STOCKHOLMFILE.close()

    return 0


if __name__ == '__main__':
    exit(main())
