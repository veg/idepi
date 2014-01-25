
from __future__ import division, print_function

from math import modf
from operator import itemgetter
from random import randint

from .._common import get_noise, sanitize_seq
from ._randomsequences import DumbRandomSequences, MarkovRandomSequences
from ._simulatedepitope import SimulatedEpitope
from ..util import is_refseq


__all__ = ['Simulation', 'DumbSimulation', 'MarkovSimulation']


def assign_class_by_percentile(seq_table, epi_def, percentile):
    vals = []
    for row in seq_table.rows():
        if is_refseq(row.id):
            continue
        vals.append(epi_def.evaluate(row.seq, get_noise(row.id)))
    proportion = calculate_percentile(vals, percentile)

    for row in seq_table.rows():
        if is_refseq(row.id):
            continue
        row.id = '|||%.3f' % epi_def.evaluate(row.seq, get_noise(row.id), proportion)

# as determined by NIST (http://www.itl.nist.gov/div898/handbook/prc/section2/prc252.htm)
def calculate_percentile(a, p):
    b = sorted(a)
    f, s = modf( p * (len(b) + 1) )
    l = int(s) - 1
    u = l + 1
    if f == 0.:
        f = 1.
    if l < 0:
        return b[0]
    elif u >= len(b):
        return b[len(b) - 1]
    return (f * b[l] + (1 - f) * b[u])

def random_column_subset(size, columns):
    col_subset = []
    assert(len(columns) > 0)
    while len(col_subset) < size:
        c_ = columns[randint(0, len(columns) - 1)]
        if c_ not in col_subset:
            col_subset.append(c_)
    return col_subset


class Simulation:
    SEQUENCE = 0
    TARGET = 1
    EPITOPE = 2
    DUMB = 3

    VALUES = (0, 1, 2, 3)


class BaseSimulation(Simulation):

    def __init__(self, mode, runs):
        self.mode = mode
        self.runs = runs

    # def simulate_epitope(self, alignment, alphabet, colnames, size, percentile, kernel_func=None):
    def simulate_epitope(self, seq_table, alphabet, column_names, size, percentile, kernel_func=None):
        if self.mode != Simulation.EPITOPE:
            return None

        if not seq_table.loaded:
            seq_table.fill_columns()

        alphdict, alphnames = alphabet.dict(), alphabet.names()
        alphabet_len = len(set(alphdict.values()))
        positions = {}

        if size > 0.2 * seq_table.num_columns:
            raise RuntimeError('we do not suggest or support simulated epitopes larger than 20% of your alignment length')

        # generate a approximately uniformly random epitope from the available nucleotides at each position (ensure that existing sequences match the epitope)
        while len(positions) < size:
            # only grab as many new positions as we need (size - len(positions))
            new_positions = random_column_subset(size - len(positions), list(seq_table.columns.keys()))
            for i in new_positions:
                # if we already have that position, delete it and skip
                if i in positions:
                    continue
                # this is some trickery to avoid ambiguous positions and spaces
                vals = [(sanitize_seq(x[0], alphdict), x[1]) for x in seq_table.columns[i].counts().items() if sanitize_seq(x[0], alphdict) not in ('X', '-')]
                # if only Xs, then ignore this position
                if len(vals) == 0:
                    continue
                count = sum([x[1] for x in vals])
                vals = sorted([(x[0], 1.0 * x[1] / count) for x in vals], key = itemgetter(1), reverse = True)
                # don't bother with this uniform bullshit, just assume the most common is the epitope
                positions[i] = vals[0][0]
                # find the value whose contribution to the cdf bounds our uniformly random value (r_) , then stop
                # r_ = random()
                # cdf = 0.
                # for j in vals:
                    # cdf += j[1]
                    # if cdf > r_:
                        # positions[i] = j[0]
                        # break

        # this formula should generate the correct position names for the epitope
        position_names = [column_names[k * alphabet_len + alphnames.index(v)] for k, v in positions.items()]

        epi_def = SimulatedEpitope(positions, position_names, alphabet, kernel_func)

        assign_class_by_percentile(seq_table, epi_def, percentile)

        return epi_def

    def generate_sequences(self, *args, **kwargs):
        raise RuntimeError('Simulation is not intended to be used directly. Please use DumbSimulation or MarkovSimulation')


class DumbSimulation(BaseSimulation):

    def __init__(self, mode, runs, refseq):
        self.__refseq = refseq
        super(DumbSimulation, self).__init__(mode, runs)

    def generate_sequences(self, N, idfmt, noise, mutation_rate, alphabet):
        return DumbRandomSequences(self.__refseq, N=N, idfmt=idfmt, noise=noise, rate=mutation_rate, alphabet=alphabet)


class MarkovSimulation(BaseSimulation):

    def __init__(self, mode, runs, refmsa):
        self.__refmsa = refmsa
        super(MarkovSimulation, self).__init__(mode, runs)

    def generate_sequences(self, N, idfmt, noise, mutation_rate, alphabet):
        return MarkovRandomSequences(self.__refmsa, N=N, idfmt=idfmt, noise=noise, rate=mutation_rate, alphabet=alphabet)
