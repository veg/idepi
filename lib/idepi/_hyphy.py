
from os import getcwd
from sys import stderr

import HyPhy as hp


__all__ = ['HyPhy']


class HyPhy(object):
    __DEFAULT_NUM_CPUS = 2

    MATRIX = hp.THYPHY_TYPE_MATRIX
    NUMBER = hp.THYPHY_TYPE_NUMBER
    STRING = hp.THYPHY_TYPE_STRING

    def __init__(self, num_cpus=None):

        if num_cpus is None:
            try:
                from multiprocessing import cpu_count
                num_cpus = cpu_count()
            except:
                num_cpus = self.__DEFAULT_NUM_CPUS

        self.__instance = hp._THyPhy(getcwd(), num_cpus)

    def execute(self, batchfile, arguments=None, quiet=True):
        execstr = 'ExecuteAFile("%s")' % batchfile if (arguments is None or len(arguments) < 1) else \
                  'ExecuteAFile("%s", { %s })' % (batchfile, ', '.join('"%d": "%s"' % (i, arguments[i]) for i in xrange(len(arguments))))
        result = self.__instance.ExecuteBF(execstr)
        out = self.__instance.GetStdout().sData
        err = self.__instance.GetErrors().sData
        war = self.__instance.GetWarnings().sData

        if not quiet:
            if war.strip() != '':
                print >> stderr, war

        if err.strip() != '':
            raise RuntimeError(err)

        return out

    def retrieve(self, variable, type):
        _res = self.__instance.AskFor(variable)
        if type not in (HyPhy.MATRIX, HyPhy.NUMBER, HyPhy.STRING):
            raise ValueError('Unknown type supplied: please use one of PhyloFilter.{MATRIX,NUMBER,STRING}')
        if (self.__instance.CanCast(_res, type)):
            res = self.__instance.CastResult(_res, type)
            if type == HyPhy.STRING:
                return res.castToString().sData
            elif type == HyPhy.NUMBER:
                return res.castToNumber().nValue
            elif type == HyPhy.MATRIX:
                return res.castToMatrix()
            else:
                # dead code, we assume
                assert(0)
        else:
            raise RuntimeError('Cast failed in HyPhy, assume an incorrect type was supplied for variable `%s\'' % variable)

