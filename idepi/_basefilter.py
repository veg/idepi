
from _alphabet import Alphabet

__all__ = ['BaseFilter']


class BaseFilter(object):

    def __init__(self):
        raise RuntimeError('You cannot use BaseFilter directly. Use one of its subclasses NaiveFilter or PhyloFilter.')

    @staticmethod
    def _labels(seqrecords, alphabet, ref_id_func, ignore_idxs=set()):
        ref = None
        for r in seqrecords:
            if apply(ref_id_func, (r.id,)):
                ref = str(r.seq)

        if ref is None:
            raise RuntimeError('No reference sequence found, aborting')

        labels = []
        c, col, ins = 0, 0, 0
        for p in ref:
            if p not in Alphabet.SPACE:
                c += 1
                ins = 0
            else:
                ins += 1
            for i in xrange(len(alphabet)):
                if col in ignore_idxs:
                    continue
                insert = base_26_to_alph(base_10_to_n(ins, BASE_ALPH))
                labels.append('%s%d%s%s' % ('' if insert != '' else p.upper(), c, insert, alphabet[i]))
                col += 1

        return labels
