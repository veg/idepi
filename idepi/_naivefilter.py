

class NaiveFilter(object):
    __DEFAULT_MIN_CONSERVATION = 1. # 100%
    __DEFAULT_MAX_CONSERVATION = 1. # 100%

    def __init__(self, seqrecords, mincons=None, maxcons=None):

        if mincons is None:
            mincons = NaiveFilter.__DEFAULT_MIN_CONSERVATION
        if maxcons is None:
            maxcons = NaiveFilter.__DEFAULT_MAX_CONSERVATION

        self.__mincons = mincons
        self.__maxcons = maxcons

        self.__seqrecords = seqrecords

    def names():
        pass
