
from __future__ import division, print_function

from collections import Iterable
from logging import getLogger
from operator import itemgetter
from sys import maxsize

from numpy import median as np_median, sum as np_sum, zeros

from ..logging import IDEPI_LOGGER


__all__ = ['Labeler']


class Labeler:

    def __init__(self, label, skip, discretize, autobalance=False):
        self.__discretize = discretize
        self.__label = label
        self.__skip = (lambda _: False) if skip is None else skip
        self.__autobalance = autobalance

    def label(self, alignment, count=None):
        return self(alignment, count)

    def __call__(self, alignment, count=None, globber=None):
        if count is None:
            count = maxsize

        if self.__autobalance:
            dtype = int
            discretize = lambda _: 1
        elif self.__discretize is None:
            dtype = float
            discretize = lambda x: x
        else:
            dtype = int
            discretize = self.__discretize

        if globber is None:
            allvals = [None] * len(alignment)
        else:
            allvals = [None] * len(globber)

        itertest = False
        label = None
        i = 0
        for seq in alignment:
            if self.__skip(seq) or i >= count:
                continue
            elif label is None:
                if isinstance(self.__label(seq), Iterable):
                    itertest = True
                    label = self.__label
                else:
                    label = lambda s: (self.__label(s),)
            if globber is None:
                r = i
            else:
                r = globber[seq.id]
            labels = sorted(label(seq), reverse=True)
            values = [discretize(x) for x in labels]
            if None in values:
                raise RuntimeError(
                    "unable to discretize label '%s'" % "','".join(labels)
                    )
            allvals[r] = list(zip(labels, values))
            i += 1

        if globber is None:
            allvals = allvals[:i]

        size = len(allvals)
        y = zeros((size,), dtype=dtype)

        median = None
        # try to balance the data
        if self.__autobalance:
            vals = allvals
            classes = old_classes = None
            median = np_median([row[0][0] for row in vals if len(row) == 1])
            ratio = old_ratio = None

            # only do this a thousand times at most
            iteration = 0
            while iteration < 1000:
                iteration += 1
                vals = [[(x, 1 if x >= median else -1) for x, c in row] for row in vals]
                classes = [set(c for _, c in row) for row in vals]
                ratio = sum(1. if 1 in s else 0. for s in classes if len(s) == 1)
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
                        max(vals[i], key=itemgetter(0))[0] if classes[i] > 0
                        else min(vals[i], key=itemgetter(0))[0]
                        for i in range(len(classes)) if len(classes[i]) == 1
                        ])

            # median doesn't change in the very last iteration of the while loop above
            discretize = lambda x: 1 if x >= median else -1

            # update allvals to only have the 1 entry for those we've already decided on,
            # the discrete value is correct because median isn't updated when we break
            # out of the loop above
            allvals = [
                [
                    max(vals[i], key=itemgetter(0)) if classes[i] > 0
                    else min(vals[i], key=itemgetter(0))
                    ] if len(classes[i]) == 1
                else vals[i]
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

        return y, median
