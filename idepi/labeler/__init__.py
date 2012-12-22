
from __future__ import division, print_function

from collections import Iterable
from logging import getLogger
from operator import itemgetter
from sys import maxsize

from numpy import median as np_median, sum as np_sum, zeros

from ..logging import IDEPI_LOGGER


__all__ = ['Labeler']


class Labeler:

    def __init__(self, label_func=None, skip_func=None, discretize_func=None, autobalance=False):
        self.__dfn = discretize_func
        self.__efn = label_func
        self.__sfn = skip_func
        self.__auto = autobalance

    def __call__(self, alignment, count=None):
        return Labeler.__label(alignment, count, self.__efn, self.__sfn, self.__dfn, self.__auto)

    def label(self, alignment, count=None):
        return Labeler.__label(alignment, count, self.__efn, self.__sfn, self.__dfn, self.__auto)

    @staticmethod
    def __label(alignment, count, label, skip, discretize, autobalance):
        if count is None:
            count = maxsize

        itertest = None
        median = None
        skipped = [False] * len(alignment)

        if autobalance:
            discretize = lambda _: True

        if discretize is None:
            dtype = float
            discretize = lambda x: x
        else:
            dtype = int

        if skip is None:
            skip = lambda _: False
        else:
            for i, row in enumerate(alignment):
                if skip(row) or i >= count:
                    skipped[i] = True
                elif itertest is None:
                    itertest = isinstance(label(alignment[i]), Iterable)

        size = len(alignment) - sum(skipped)
        y = zeros((size,), dtype=dtype)

        if size == 0:
            raise RuntimeError('skipping EVERYTHING in Labeler?!')

        allvals = None
        if itertest:
            allvals = [
                [
                    (x, discretize(x)) for x in sorted(label(row), reverse=True)
                ] for i, row in enumerate(alignment) if not skipped[i]
            ]
        else:
            allvals = [
                [
                    (x, discretize(x)) for x in [label(row)]
                ] for i, row in enumerate(alignment) if not skipped[i]
            ]

        # try to balance the data
        if autobalance:
            vals = allvals
            classes = old_classes = None
            median = np_median([row[0][0] for row in vals if len(row) == 1])
            ratio = old_ratio = None

            # only do this ten thousand times at most
            iteration = 0
            while iteration < 1000:
                iteration += 1
                vals = [[(x, x >= median) for x, c in row] for row in vals]
                classes = [set(c for _, c in row) for row in vals]
                ratio = sum(1. if True in s else 0. for s in classes if len(s) == 1)
                # break if we stop changing class evaluation
                if old_classes is not None and classes == old_classes:
                    break
                # if the old ratio resulted in a better partitioning,
                # then keep that result and terminate, otherwise keep trying.
                elif old_ratio is not None and abs(0.5 - old_ratio) < abs(0.5 - ratio):
                    classes = old_classes
                    ratio = old_ratio
                    break
                else:
                    old_classes = classes
                    # set contains at least one and at most two values: True or False
                    # if True is not in set, then False is, thus the median is computed:
                    # take the maximum value if len(set) == 1 and True in set otherwise
                    # take the minimum value if len(set) == 1 and True not in set otherwise
                    # don't include it
                    median = np_median([
                        max(vals[i], key=itemgetter(0))[0] if True in classes[i] else
                        min(vals[i], key=itemgetter(0))[0]
                        for i in range(len(classes)) if len(classes[i]) == 1
                    ])

            # median doesn't change in the very last iteration of the while loop above
            discretize = lambda x: x >= median

            # update allvals to only have the 1 entry for those we've already decided on,
            # the discrete value is correct because median isn't updated when we break
            # out of the loop above
            allvals = [
                [
                    max(vals[i], key=itemgetter(0)) if True in classes[i] else
                    min(vals[i], key=itemgetter(0))
                ] if len(classes[i]) == 1 else vals[i]
                for i in range(len(classes))
            ]

        if itertest:
            ambigs = {}
            for i, values in enumerate(allvals):
                classes = set(v for k, v in values)
                if len(classes) > 1:
                    ambigs[i] = values
                else:
                    y[i] = classes.pop()

            classavg = np_sum(y) / (size - len(ambigs))
            pos = min(max(int((0.5 - classavg) * size + 0.5), 0), len(ambigs))

            log = getLogger(IDEPI_LOGGER)
            log.debug('found %d ambiguous in %d records%s' % (len(ambigs), size,
                    ', %d to be interpreted positively' % pos if pos > 0 else ''
            ))

            for i in range(pos):
                # kv is key-value,
                # so kv[1] is the revsorted list [(ic50, klass), ...],
                # and kv[1][0][0] is the largest ic50 value for key k
                kv = max(ambigs.items(), key=lambda kv: kv[1][0][0])
                idx, klass = kv[0], kv[1][0][1]
                y[idx] = klass
                del ambigs[idx]
            # remaining are to be left at 0
        else:
            for i, values in enumerate(allvals):
                y[i] = values[0][1] # 0 automagically exists, and 1 refers to the discretized value

        # set the values between -1 and 1
        y *= 2
        y -= 1

        return y, median if autobalance else None
