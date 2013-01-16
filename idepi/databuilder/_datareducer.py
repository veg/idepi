
from numpy import hstack


__all__ = ['DataReducer']


class DataReducer:

    def __init__(self, *args):
        self.__builders = args
        self.__labels = []
        self.__length = 0
        for builder in self.__builders:
            self.__labels.extend(builder.labels)
            self.__length += len(builder)

    def __call__(self, *args, **kwargs):
        return hstack(
            tuple(
                builder(*args, **kwargs)
                for builder in self.__builders
                )
            )

    def __len__(self):
        return self.__length

    @property
    def labels(self):
        return self.__labels
