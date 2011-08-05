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

import numpy as np

from _normalvalue import NormalValue


__all__ = ['PerfStats']


class PerfStats(object):
    CONTINUOUS  = 0
    DISCRETE    = 1

    ACCURACY            = 0
    PPV, PRECISION      = 1, 1
    NPV                 = 2
    SENSITIVITY, RECALL = 3, 3
    SPECIFICITY, TNR    = 4, 4
    FSCORE              = 5
    MINSTAT             = 6
    
    RBARSQUARED         = 7
    RSQUARED            = 8
    RMSE                = 9

    def __init__(self, mode=None):
        if mode is None:
            mode = PerfStats.DISCRETE
        
        if mode not in (PerfStats.CONTINUOUS, PerfStats.DISCRETE):
            raise RuntimeError('PerfStat mode must be one of either PerfStats.CONTINUOUS or PerfStats.DISCRETE')

        self.mode = mode

        if self.mode == PerfStats.DISCRETE:
            self.accuracy = NormalValue(float)
            self.ppv = NormalValue(float)
            self.npv = NormalValue(float)
            self.sensitivity = NormalValue(float)
            self.specificity = NormalValue(float)
            self.fscore = NormalValue(float)
            self.minstat = NormalValue(float)
        else:
            self.rbar2 = NormalValue(float)
            self.r2 = NormalValue(float)
            self.rmse = NormalValue(float)

    def append(self, truth, preds, weights=None):
        if self.mode == PerfStats.CONTINUOUS:
            if weights is None:
                raise ValueError('weights must not be None for PerfStats.CONTINUOUS')
            rbar2, r2, rmse = PerfStats.calcstat_continuous(truth, preds, weights)
            self.rbar2.append(rbar2)
            self.r2.append(r2)
            self.rmse.append(rmse)
        else:
            acc, ppv, npv, sen, spe, fsc, mst = PerfStats.calcstat_discrete(truth, preds)
            self.accuracy.append(acc)
            self.ppv.append(ppv)
            self.npv.append(npv)
            self.sensitivity.append(sen)
            self.specificity.append(spe)
            self.fscore.append(fsc)
            self.minstat.append(mst)

    @staticmethod
    def calcstat_continuous(y, yhat, w):
        nperr = np.seterr(divide='ignore')
        sse = sum(pow(y - yhat, 2.0))
        ybar = np.mean(y)
        sst = sum(pow(y - ybar, 2.0))
        r2 = 1.0 - (sse / sst)
        nless1 = len(y) - 1
        p = len([1 for i in w if i != 0.0]) - 1 # - 1 to avoid counting the constant term
        mse = sse / (nless1 - p) # count the the full N
        rmse = np.sqrt(mse)
        rbar2 = 1.0 - (1.0 - r2) * nless1 / (nless1 - p)

        np.seterr(**nperr)

        return rbar2, r2, rmse

    @staticmethod
    def calcstat_discrete(truth, preds):
        tp, tn, fp, fn = PerfStats.ystoconfusionmatrix(truth, preds)

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

        return acc, ppv, npv, sen, spe, fsc, mst

    def get(self, stat):
        if self.mode == PerfStats.CONTINUOUS:
            if stat == PerfStats.RMSD:
                return self.rmsd
            else:
                pass
        else: # PerfStats.DISCRETE
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
                pass
        raise ValueError('No such statistic exists here.')

    def tolist(self):
        if self.mode == PerfStats.CONTINUOUS:
            return [
                (u'R\u0304\u00b2', self.rbar2),
                (u'R\u00b2', self.r2),
                (u'RMSE', self.rmse)
            ]
        else: # PerfStats.DISCRETE
            return [
                (u'Accuracy', self.accuracy),
                (u'PPV', self.ppv),
                (u'NPV', self.npv),
                (u'Sensitivity', self.sensitivity),
                (u'Specificity', self.specificity),
                (u'F-score', self.fscore),
                (u'Minstat', self.minstat)
            ]

    def todict(self):
        return dict(PerfStats.tolist(self))

    @staticmethod
    def ystoconfusionmatrix(truth, preds):
        tps = truth > 0.
        tns = truth <= 0.
        pps = preds > 0.
        pns = preds <= 0.
                                                               # true pos    true neg    false pos   false neg
        tp, tn, fp, fn = map(lambda a: np.sum(np.multiply(*a)), [(tps, pps), (tns, pns), (tns, pps), (tps, pns)])

        return tp, tn, fp, fn

    def __str__(self):
        return str(PerfStats.todict(self))
