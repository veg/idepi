#!/usr/bin/env python2.7

from os.path import basename, exists
from sys import argv as sys_argv, exit as sys_exit, stderr

from Bio import SeqIO

from idepi import Alphabet, Ancestors, id_to_float, is_HXB2


def main(argv=sys_argv):
    name = basename(argv.pop(0))

    try:
        features = set(f.upper() for f in argv.pop().split(','))
        assert(len(argv) == 1 and exists(argv[0]) and len(features))
    except:
        print >> stderr, 'usage: %s ALIGNMENT FEATURES' % name
        sys_exit(-1)

    with open(argv[0]) as fh:
        seqrecords = [r for r in SeqIO.parse(fh, 'stockholm')]

    alph = Alphabet(Alphabet.AMINO)

    refseq = None
    for i, r in enumerate(seqrecords):
        if is_HXB2(r.id):
            refseq = str(seqrecords[i].seq).upper()
            del seqrecords[i]
            break

    if refseq is None:
        raise RuntimeError('no reference sequence was found in the alignment, aborting')

    colnames = []
    colnum = 1
    for i in xrange(len(refseq)):
        if refseq[i] not in '._-':
            colname = refseq[i] + str(colnum)
            colnum += 1
            if colname in features:
                colnames.append((i, colname)) 

    tree, alignment = Ancestors()(seqrecords)

    for r in alignment:
        # labels has length of colnames plus the ic50
        labels = [None] * (len(colnames) + 1)
        i = 1
        for idx, colname in colnames:
            labels[i] = colname + r.seq[idx]
            i += 1
        try:
            labels[0] = '%.3g' % id_to_float(r.id)
        except ValueError:
            labels[0] = ''
        tree = tree.replace(r.id, '_'.join(labels))

    print tree

    return 0

if __name__ == '__main__':
    sys_exit(main())
