
from os import getcwd
from sys import stderr

import HyPhy as hp


__all__ = ['HyPhy']


class HyPhy(object):
    __DEFAULT_NUM_CPUS = 2

    def __init__(self, num_cpus=None):

        if num_cpus is None:
            try:
                from multiprocessing import cpu_count
                num_cpus = cpu_count()
            except:
                num_cpus = self.__DEFAULT_NUM_CPUS

        self.__instance = hp._THyPhy(getcwd(), num_cpus)

    def execute(self, batchfile, quiet=False):
        result = self.__instance.ExecuteBF('ExecuteAFile("%s")' % batchfile)
        out = self.__instance.GetStdout()
        err = self.__instance.GetErrors()
        war = self.__instance.GetWarnings()

        if not quiet:
            if war.strip() != '':
                print >> stderr, war

        if err.strip() != '':
            raise RuntimeError(err)

        return out 
