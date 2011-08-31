
from sys import maxint

from numpy import zeros


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
            count = maxint

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
                if apply(skip, (row.id,)) or i > count:
                    skipped += 1

        i = 0
        y = zeros((len(alignment) - skipped,), dtype=dtype)

        off = 0
        for i, row in enumerate(alignment):
            if apply(skip, (row.id,)):
                off += 1
                pass
            else:
                y[i - off] = apply(discretize, (apply(extract, (row.id,)),))
            if i > count:
                break

        return y
