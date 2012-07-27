#!/usr/bin/env python

from __future__ import division, print_function

from os.path import basename, exists
from re import compile as re_compile, sub as re_sub
from sys import argv as sys_argv, exit as sys_exit, stderr

from numpy import mean

from idepi import PhyloGzFile, seqrecord_get_ic50s


def main(argv=sys_argv):
    name = basename(argv.pop(0))
    numeric = re_compile(r'[^0-9]+')

    try:
        features = set(int(numeric.sub('', f)) for f in argv.pop().split(','))
        assert(len(argv) == 1 and exists(argv[0]) and len(features))
    except:
        print('usage: %s TREE FEATURES' % name, file=stderr)
        sys_exit(-1)

    tree, alignment, colnames, _ = PhyloGzFile.read(argv[0])

    icolnames = [(idx, colname) for idx, colname in enumerate(colnames) if int(numeric.sub('', colname)) in features]

    for r in alignment:
        # labels has length of icolnames plus the ic50
        labels = [None] * (len(icolnames) + 1)
        i = 1
        for idx, colname in icolnames:
            if len(colnames) > 1:
                labels[i] = colname + r.seq[idx]
            else:
                labels[i] = r.seq[idx]
            i += 1
        try:
            labels[0] = '%.3g' % mean(seqrecord_get_ic50s(r))
        except ValueError:
            if not (len(r.id) > 4 and r.id[:4].lower() == 'node'):
                print(r)
            labels.pop(0)
        # include the ':' here to make sure we grab the end of the label
        tree = re_sub(r'([,()])' + r.id + r'(?:_[0-9]+)?:', r'\g<1>' + '_'.join(labels) + ':', tree)

    print(tree)

    return 0

if __name__ == '__main__':
    sys_exit(main())
