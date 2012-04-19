
from collections import defaultdict
from operator import itemgetter


class Counter(defaultdict):

    def __init__(self, arg):
        counts = {}
        for a in arg:
            if a not in counts:
                counts[a] = 0
            counts[a] += 1
        super(Counter, self).__init__(lambda: 0, counts)

    def elements(self):
        lst = []
        for a, c in self.items():
            if c > 0:
                lst.extend([a] * c)
        return lst

    def most_common(self, n=None):
        lst = sorted(self.items(), key=itemgetter(1), reverse=True)
        if n is not None and n < len(lst):
            return lst[:n]
        else:
            return lst

    def subtract(self, iterable=[]):
        s = Counter(iterable)
        for a, c in s.items():
            self[a] -= c

    @staticmethod
    def fromkeys():
        raise NotImplementedError

    def update(self, iterable):
        s = Counter(iterable)
        for a, c in s.items():
            self[a] += c
