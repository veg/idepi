#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# mrmr.py :: computes Maximum Relevance (MaxRel) and minimum
# Redundancy Maximum Relevance (mRMR) for a dataset
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

import os, sys, time

import multiprocessing as mp
import numpy as np

import platform
if platform.system().lower() == 'darwin':
    __subpath = os.path.join('python%d%d' % sys.version_info[:2], 'lib', 'python')
    # numpy must go first
    for module in ('numpy-1.5.1', 'biopython-1.56', 'cvxopt-1.1.3', 'cython-0.14.1', 'mlpy-2.2.2', 'pil-1.1.7'):
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'contrib', module, __subpath)
        if os.path.exists(path):
            sys.path.insert(0, path)

from idepi import MAXREL, MID, MIQ, Mrmr, UiThread


THRESHOLD = 0.8


def __usage():
    print >> sys.stderr, 'usage: %s MRMRFILE' % __name
    return -1


def main(argv, ui=None):

    num_features = 10

    i = 0
    while i < len(argv):
        if argv[i][0] == '-':
            opt = argv[i][1]
            if opt is 'n':
                num_features = int(argv[i+1])
                i += 2
                continue
            elif opt is '-':
                if argv[i][2:] == 'gpl':
                    print gplv2
                    return 0
            return usage()
        else:
            argv = argv[i:]
            break

    if len(argv) != 1 or not os.path.exists(argv[0]):
        return usage()

    fh = open(argv[0], 'r')
    lines = [l.strip() for l in fh]
    fh.close()

    names = lines[0].split(',')[1:]
    arr = [[int(x) for x in l.split(',')] for l in lines[1:] if l != '']
    m = np.matrix(arr, dtype=bool)

    targets = m[:, 0]
    vars = m[:, 1:]

    selector = Mrmr()

    # hax to get at 'private' method 
    maxrel, mrmr = selector._Mrmr__mrmr_selection(num_features, MID, vars, targets, threshold=THRESHOLD, ui=ui)

    print 'I(X, Y) / H(X, Y)'
    for idx, value in maxrel: 
        print '   %4d   % 5s   %6.4f' % (idx + 1, names[idx], value)
   
    print '\nI(X, Y) / H(X, Y) - I(X, X) / H(X, X) (related > %.3g)' % THRESHOLD    
    for idx, value, related in mrmr:
        print '   %4d   % 5s   %6.4f   (%s)' % (idx + 1, names[idx], value, ', '.join([(names[i], v) for i, v in related]))
    
    return 0 


if __name__ == '__main__':
    __name = os.path.basename(sys.argv[0])
    usage = __usage
    ui = UiThread()
    sys.exit(main(sys.argv[1:], ui))
