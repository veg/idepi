
from functools import partial


__all__ = ['naivefilter']


def __naivefilter(maxcons, mincons, maxgap, pos):
    # return the number of binary columns per site
    if (pos is None or
        pos.maxcount / pos.total > maxcons or
        pos.mincount / pos.total > mincons or
        pos.gapcount / (pos.total + pos.gapcount) > maxgap):
        return []
    # otherwise return true
    return sorted(pos.counts.keys())

def naivefilter(maxcons, mincons, maxgap):
    return partial(__naivefilter, maxcons, mincons, maxgap)
