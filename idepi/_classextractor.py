
from numpy import zeros


__all__ = ['ClassExtractor']


class ClassExtractor(object):

    def __init__(self, extract_func=None, skip_func=None, discretize_func=None):
        self.__dfn = discretize_func
        self.__efn = extract_func
        self.__sfn = skip_func

    def extract(self, alignment):
        return ClassExtractor.__extract(alignment, self.__efn, self.__sfn, self.__dfn)

    @staticmethod
    def __extract(alignment, extract, skip, discretize):
        if discretize is None:
            dtype = float
            discretize = lambda x: x
        else:
            dtype = bool

        skipped = 0
        if skip is None:
            skip = lambda _: False
        else:
            for row in alignment:
                if apply(skip, (row.id,)):
                    skipped += 1

        i = 0
        y = zeros((len(alignment) - skipped,), dtype=dtype)

        for row in alignment:
            if apply(skip, (row.id,)):
                continue
            else:
                y[i] = apply(discretize, (apply(extract, (row.id,)),))
                i += 1
        return y
