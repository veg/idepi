
from collections import Iterable
from sys import maxsize

from numpy import sum as np_sum, zeros


__all__ = ['ClassExtractor']


class ClassExtractor(object):

    def __init__(self, extract_func=None, skip_func=None, discretize_func=None):
        self.__dfn = discretize_func
        self.__efn = extract_func
        self.__sfn = skip_func

    def extract(self, alignment, count=None):
        return ClassExtractor.__extract(alignment, count, self.__efn, self.__sfn, self.__dfn)

    @staticmethod
    def __extract(alignment, count, extract, skip, discretize):
        if count is None:
            count = maxsize

        skipped = 0
        if discretize is None:
            dtype = float
            discretize = lambda x: x
        else:
            dtype = bool

        if skip is None:
            skip = lambda _: False
        else:
            for i, row in enumerate(alignment):
                if skip(row) or i >= count:
                    skipped += 1

        size = len(alignment) - skipped
        y = zeros((size,), dtype=dtype)

        # try to balance the data
        if isinstance(extract(alignment[0]), Iterable):
            ambigs = {}
            off = 0
            for i, row in enumerate(alignment):
                if skip(row) or i >= count:
                    off += 1
                else:
                    values = [(x, discretize(x)) for x in sorted(extract(row), reverse=True)]
                    classes = set(v for k, v in values)
                    if len(classes) > 1:
                        ambigs[i] = off, values
                    else:
                        y[i - off] = classes.pop()
                if i >= count:
                    break
            classavg = np_sum(y) / (size - len(ambigs))
            pos = max(int((0.5 - classavg) * len(ambigs) + 0.5), 0)
            for i in range(min(pos, len(ambigs))):
                # kv is key-value, so kv[1] is value,
                # kv[1][1] is the revsorted list [(ic50, klass), ...],
                # and kv[1][1][0][0] is the largest ic50 value for key k
                kv = max(ambigs.items(), key=lambda kv: kv[1][1][0][0])
                idx, off, klass = kv[0], kv[1][0], kv[1][1][0][1]
                y[idx - off] = klass
                del ambigs[idx]
            # remaining are to be left at 0
        else:
            off = 0
            for i, row in enumerate(alignment):
                if skip(row) or i >= count:
                    off += 1
                else:
                    y[i - off] = discretize(extract(row))

        return y
