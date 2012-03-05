
from logging import getLogger
from time import clock

from numpy import zeros

from ._alphabet import Alphabet
from ._logging import IDEPI_LOGGER
from ._util import BASE_ALPH, base_10_to_n, base_26_to_alph


class BaseFilter(object):

    def __init__(self):
        raise RuntimeError('You cannot use BaseFilter directly. Use one of its subclasses NaiveFilter or PhyloFilter.')

    # refseq_offs has for keys 0-indexed positions into the trimmed refseq
    # and for values the number of trimmed positions occuring immediately
    # before the key-index. This so the colfilter can properly name the
    # positions according to the full reference sequence
    @staticmethod
    def _colnames(alignment, alphabet, ref_id_func, refseq_offs, ignore_idxs):
        ref = None
        for r in alignment:
            if ref_id_func(r):
                ref = str(r.seq)

        if ref is None:
            raise RuntimeError('No reference sequence found, aborting')

        b = clock()
        colnames = []
        c, col, ins = 0, 0, 0
        for i, p in enumerate(ref):
            if p not in Alphabet.SPACE:
                c += 1
                ins = 0
            else:
                ins += 1
            if i in refseq_offs:
                # if we're in refseq_offs, then skip
                c += refseq_offs[i]
            for j in range(len(alphabet)):
                if col in ignore_idxs:
                    pass
                else:
                    insert = base_26_to_alph(base_10_to_n(ins, BASE_ALPH))
                    colnames.append('%s%d%s%s' % ('' if insert != '' else p.upper(), c, insert, alphabet[j]))
                col += 1

        getLogger(IDEPI_LOGGER).debug('finished labeling columns, took %.3f' % (clock() - b))

        return colnames
