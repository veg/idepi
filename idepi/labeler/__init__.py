
from __future__ import division, print_function

from json import loads as json_loads

from numpy import zeros


__all__ = [
    'Labeler',
    'expression',
    'skipper',
    ]


def expression(label, seqrecord):
    desc = json_loads(seqrecord.description)
    try:
        return eval(label, {}, desc['values'])
    except NameError:
        return None


def skipper(is_refseq, subtypes, seqrecord):
    if is_refseq(seqrecord):
        return True
    if not subtypes:
        return False
    desc = json_loads(seqrecord.description)
    try:
        return desc['subtype'] not in subtypes
    except KeyError:
        return True


class Labeler:

    def __init__(self, label, skip):
        self.__label = label
        self.__skip = (lambda _: False) if skip is None else skip

    def label(self, alignment):
        return self(alignment)

    def __call__(self, alignment):

        size = len(alignment)

        labels = [None for _ in range(size)]

        alignment_ = alignment[:0]
        i = 0
        for seq in alignment:
            if self.__skip(seq):
                continue
            label = self.__label(seq)
            # skip if we have no label at all
            if label is None:
                continue
            labels[i] = label
            alignment_.append(seq)
            i += 1

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
