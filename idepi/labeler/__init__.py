
from __future__ import division, print_function

from json import loads as json_loads

from numpy import median, zeros

from Bio.Align import MultipleSeqAlignment


__all__ = [
    'Labeler',
    'expression'
    ]



def expression(seqrecord, label):
    desc = json_loads(seqrecord.description)
    try:
        return eval(label, {}, desc['values'])
    except NameError:
        return None


class Labeler:

    def __init__(self, label, skip, autobalance=False):
        self.__label = label
        self.__skip = (lambda _: False) if skip is None else skip
        self.__autobalance = autobalance

    def label(self, alignment, globber=None):
        return self(alignment, globber)

    def __call__(self, alignment, globber=None):

        size = len(alignment) if globber is None else len(globber)

        labels = [None for _ in range(size)]

        alignment_ = MultipleSeqAlignment(
            [],
            alphabet=alignment._alphabet
            )
        i = 0
        for seq in alignment:
            if self.__skip(seq):
                alignment_.append(seq)
                continue
            label = self.__label(seq)
            # skip if we have no label at all
            if label is None:
                continue
            if globber is None:
                r = i
            else:
                r, _ = globber[seq.id]
            labels[r] = label
            alignment_.append(seq)
            i += 1

        if globber is None:
            size = i
            labels = labels[:size]

        dtype = int
        for label in labels:
            type_ = type(label)
            # bools are ints
            if type_ is bool:
                type_ = int
            # raise errors is something is wonky
            if type_ not in (int, float):
                raise ValueError(
                    "unmanageable label type '{0}'".format(
                        type_.__name__
                        )
                    )
            # if even a single value is a float,
            # they're all floats
            if dtype is int and type_ is float:
                dtype = float

        y = zeros((size,), dtype=dtype)

        # try to balance the data using median
        median_ = None
        if self.__autobalance:
            median_ = median(labels)
            labels = [label > median_ for label in labels]

        i = 0
        for label in labels:
            if label is None:
                continue
            # convert to +1/-1 if we're bool
            if type(label) is bool:
                y[i] = 1 if label else -1
            else:
                y[i] = label
            i += 1

        return alignment_, y, median_
