
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
        if discretize is None:
            dtype = float
            discretize = lambda x: x
        else:
            dtype = bool

        skipped = 0
        if skip is None:
            skip = lambda _: False
        else:
            for i, row in enumerate(alignment):
                if apply(skip, (row.id,)) or i > count:
                    skipped += 1

        i = 0
        y = zeros((len(alignment) - skipped,), dtype=dtype)

        for j, row in enumerate(alignment):
            if apply(skip, (row.id,)):
                pass
            else:
                y[i] = apply(discretize, (apply(extract, (row.id,)),))
                i += 1
            if j > count:
                break

        # this should be true
        assert(j) == len(y)

        return y
