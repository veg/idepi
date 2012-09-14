
from numpy import hstack


__all__ = ['DataReducer']


class DataReducer(object):

    def __init__(self, *args):
        self.__dbs = args
        self.__labels = []
        self.__length = 0
        for db in self.__dbs:
            self.__labels.extend(db.labels)
            self.__length += len(db)

    def __call__(self, alignment, refidx=None):
        return hstack(tuple(db(alignment, refidx) for db in self.__dbs))

    def __len__(self):
        return self.__length

    @property
    def labels(self):
        return self.__labels
