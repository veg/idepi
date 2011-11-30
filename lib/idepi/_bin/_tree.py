#!/usr/bin/env python2.7

from gzip import GzipFile
from os.path import basename, exists
from re import compile as re_compile
from sys import argv as sys_argv, exit as sys_exit, stderr

try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO

from Bio import SeqIO

from idepi import id_to_float


def main(argv=sys_argv):
    name = basename(argv.pop(0))
    numeric = re_compile(r'[^0-9]+')

    try:
        features = set(int(numeric.sub('', f)) for f in argv.pop().split(','))
        assert(len(argv) == 1 and exists(argv[0]) and len(features))
    except:
        print >> stderr, 'usage: %s TREE FEATURES' % name
        sys_exit(-1)

    tree = None
    alignment = None
    refseq = None
    buf = None
    with GzipFile(argv[0]) as fh:
        for line in (l.strip() for l in fh):
            sep = '' if len(line) < 5 else line[:5].lower()
            if sep.find('end') == 0 and buf is not None:
                typ = line.split()[1].lower()
                if typ.find('newick') >= 0:
                    tree = buf.getvalue().strip()
                elif typ.find('refseq') >= 0:
                    refseq = buf.getvalue().strip()
                elif typ.find('fasta') >= 0:
                    buf.seek(0)
                    alignment = [r for r in SeqIO.parse(buf, 'fasta')]
                buf.close()
                buf = None
            elif sep.find('begin') == 0 and buf is None:
                buf = StringIO()
            elif buf is not None:
                buf.write(line)
                buf.write('\n')

    if buf is not None:
        buf.close()

    colnames = [None] * len(features)
    colnum = 1
    j = 0
    for i in xrange(len(refseq)):
        if refseq[i] not in '._-':
            colname = refseq[i] + str(colnum)
            if colnum in features:
                colnames[j] = (i, colname)
                j += 1
            colnum += 1

    for r in alignment:
        # labels has length of colnames plus the ic50
        labels = [None] * (len(colnames) + 1)
        i = 1
        for idx, colname in colnames:
            if len(colnames) > 1:
                labels[i] = colname + r.seq[idx]
            else:
                labels[i] = r.seq[idx]
            i += 1
        try:
            labels[0] = '%.3g' % id_to_float(r.id)
        except ValueError:
            labels.pop(0)
        # include the ':' here to make sure we grab the end of the label
        tree = tree.replace(r.id + ':', '_'.join(labels) + ':')

    print tree

    return 0

if __name__ == '__main__':
    sys_exit(main())
