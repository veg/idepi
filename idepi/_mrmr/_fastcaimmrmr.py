
from _fastcaim import FastCaim
from _discretemrmr import DiscreteMrmr


__all__ = ['FastCaimMrmr']


def FastCaimMrmr(DiscreteMrmr):

    def __init__(self, *args, **kwargs):
        self.__selected = False
        self.__fc = FastCaim()
        super(FastCaimMrmr, self).__init__(*args, **kwargs)

    def select(self, x, y):
        self.__selected = False
        self.__fc.learn(x, y)
        x = self.__fc.discretize(x)
        r = DiscreteMrmr.select(self, x, y)
        self.__selected = True
        return r

    def subset(self, x):
        if not self.__selected:
            raise RuntimeError('No FastCaimMrmr model computed.')
        x = self.__fc.discretize(x)
        return DiscreteMrmr.subset(self, x)
