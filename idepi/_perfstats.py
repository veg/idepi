#
# idepi :: (IDentify EPItope) python libraries containing some useful machine
# learning interfaces for regression and discrete analysis (including
# cross-validation, grid-search, and maximum-relevance/mRMR feature selection)
# and utilities to help identify neutralizing antibody epitopes via machine
# learning.
#
# Copyright (C) 2011 N Lance Hepler <nlhepler@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from _normalvalue import NormalValue


__all__ = ['PerfStats']


class PerfStats(object):

    ACCURACY            = 0
    PPV, PRECISION      = 1, 1
    NPV                 = 2
    SENSITIVITY, RECALL = 3, 3
    SPECIFICITY, TNR    = 4, 4
    FSCORE              = 5
    MINSTAT             = 6

    def __init__(self):
        self.accuracy = NormalValue(float)
        self.ppv = NormalValue(float)
        self.npv = NormalValue(float)
        self.sensitivity = NormalValue(float)
        self.specificity = NormalValue(float)
        self.fscore = NormalValue(float)
        self.minstat = NormalValue(float)

    def append(self, tp, tn, fp, fn):
        # convert to ints otherwise numpyshit may happen and screw crap up
        tp, tn, fp, fn = int(tp), int(tn), int(fp), int(fn)

        tot = tp + tn + fp + fn
        acc = 0. if tot == 0 else float(tp + tn) / tot
        ppv = 0. if tp + fp == 0 else float(tp) / (tp + fp)
        npv = 0. if tn + fn == 0 else float(tn) / (tn + fn)
        sen = 0. if tp + fn == 0 else float(tp) / (tp + fn)
        spe = 0. if tn + fp == 0 else float(tn) / (tn + fp)
        fsc = 0. if ppv + sen == 0. else 2. * ppv * sen / (ppv + sen)
        mst = min(ppv, npv, sen, spe) # fscore can't me min as a harmonic mean

        self.accuracy.append(acc)
        self.ppv.append(ppv)
        self.npv.append(npv)
        self.sensitivity.append(sen)
        self.specificity.append(spe)
        self.fscore.append(fsc)
        self.minstat.append(mst)

    def get(self, stat):
        if stat == PerfStats.ACCURACY:
            return self.accuracy
        elif stat == PerfStats.PPV:
            return self.ppv
        elif stat == PerfStats.NPV:
            return self.npv
        elif stat == PerfStats.SENSITIVITY:
            return self.sensitivity
        elif stat == PerfStats.SPECIFICITY:
            return self.specificity
        elif stat == PerfStats.FSCORE:
            return self.fscore
        elif stat == PerfStats.MINSTAT:
            return self.minstat
        else:
            raise ValueError('No such statistic exists here.')

    def tolist(self):
        return [
            ('Accuracy', self.accuracy),
            ('PPV', self.ppv),
            ('NPV', self.npv),
            ('Sensitivity', self.sensitivity),
            ('Specificity', self.specificity),
            ('F-score', self.fscore),
            ('Minstat', self.minstat)
        ]

    def todict(self):
        return dict(PerfStats.tolist(self))

    def __str__(self):
        return str(PerfStats.todict(self))
