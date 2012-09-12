
from __future__ import division, print_function

from logging import getLogger
from time import clock

from BioExt import Counter

from ..alphabet import Alphabet
from ..logging import IDEPI_LOGGER


def consgap_skip_columns(alignment, mincons, maxcons, maxgap, refidx=None):

    skip_colidxs = set()

    # if these variables are all set to 1, then
    # there is no work to do
    if mincons < 1. or maxcons < 1. or maxgap < 1.:

        aln_len = alignment.get_alignment_length()
        b = clock()

        for i in range(aln_len):
            counts = Counter(alignment[j, i] for j in range(len(alignment)) if j != refidx)
            col = [v for k, v in counts.items() if k != Alphabet.SPACE]
            colsum_nogaps = sum(col)
            colsum_withgaps = colsum_nogaps + (counts[Alphabet.SPACE] if Alphabet.SPACE in counts else 0)
            # kill perfectly conserved or empty columns
            if (colsum_nogaps == 0. or
                min(col) / colsum_nogaps > mincons or
                max(col) / colsum_nogaps > maxcons or
                max(col) / colsum_nogaps >= 1. or
                (Alphabet.SPACE in counts and counts[Alphabet.SPACE] / colsum_withgaps > maxgap)):
                skip_colidxs.add(i)

        getLogger(IDEPI_LOGGER).debug('finished conservation and gap filtration, took %.3f' % (clock() - b))

    return skip_colidxs
