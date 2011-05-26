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

from exceptions import StandardError
from re import match
from operator import itemgetter
from os import remove
from os.path import exists
from shutil import copyfile
from subprocess import PIPE, Popen
from sys import exit, stderr, stdout
from tempfile import mkstemp

import numpy as np

from _fakepool import FakePool


__all__ = ['Mrmr', 'MAXREL', 'MID', 'MIQ']

MAXREL = 0
MID = 1
MIQ = 2


class Mrmr(object):
    __DEFAULT_THRESHOLD = 0.8

    def __init__(self, num_features=10, method=MID, threshold=None):
        if threshold is None:
            threshold = self.__DEFAULT_THRESHOLD

        self.__computed = False
        self.__colsize = 0
        self.__maxrel = None
        self.__mrmr = None
        self.method = method
        self.num_features = num_features
        self.threshold = threshold

    @staticmethod
    def __compute_mi_inner(x, vars_v, targets_t, p=None):
        p_t = float(np.sum(targets_t)) / x # p(X == t)
        p_v = np.sum(vars_v, axis=0).astype(float) / x # p(Y == v)
        p_tv = np.sum(np.multiply(targets_t, vars_v), axis=0).astype(float) / x # p(X == t, Y == v)
        mi = np.nan_to_num(np.multiply(p_tv, np.log2(p_tv / (p_t * p_v))))
        h = -np.nan_to_num(np.multiply(p_tv, np.log2(p_tv)))
    
        if p is not None:
            p.value += 1
        
        return mi, h

    @classmethod
    def __compute_mi(cls, x, vars, targets, progress=None):
    
        vcache = { True: vars == True, False: vars == False }
        tcache = { True: targets == True, False: targets == False }
        
        res = {}
        
        pool = FakePool() # mp.Pool(mp.cpu_count())
    
        for t in (True, False):
            for v in (True, False):
                res[(t, v)] = pool.apply_async(Mrmr.__compute_mi_inner, (x, vcache[v], tcache[t], progress))
        
        pool.close()
        pool.join()
        
        y = vars.shape[1]
        mi, h = np.zeros(y, dtype=float), np.zeros(y, dtype=float)
    
        for r in res.values():
            __mi, __h = r.get()
            mi = np.add(mi, __mi)
            h = np.add(h, __h)
            if progress is not None:
                progress.value += 1
    
        return mi, h

    @staticmethod
    def __compute_mi_xbar(mi, a):
        return (np.sum(mi[a, :]) - mi[a, a]) / (mi.shape[0] - 1)    
    
    @staticmethod
    def __compute_mibar(mi):
        n = mi.shape[0]
        m = n - 1
        mibar = 0.
        for x in xrange(m-1):
            for y in xrange(x+1, n-1):
                mibar += mi[x, y]
        return 2. * mibar / (m * n) 
    
    @classmethod
    def __compute_apc(cls, mi, a, b, mibar=None):
        if mibar is None:
            mibar = Mrmr.__compute_mibar(mi)
        return Mrmr.__compute_mi_xbar(mi, a) * Mrmr.__compute_mi_xbar(mi, b) / mibar

    # taken from Dunn et al 2007, 'Mutual information without the influence
    # of phylogeny or entropy dramatically improves residue contact prediction',
    # Bioinformatics (2008) 24 (3): 333-340
    @classmethod
    def __compute_mip(cls, mi, a, b, mibar=None):
        if mibar is None:
            mibar = Mrmr.__compute_mibar(mi)
        return mi[a, b] - Mrmr.__compute_apc(mi, a, b, mibar)

    @classmethod
    def __mrmr_selection(cls, num_features, method, vars, targets, threshold=None, ui=None):
        if method not in (MAXREL, MID, MIQ):
            raise ValueError('method must be one of Mrmr.MAXREL, Mrmr.MID, or Mrmr.MIQ')

        if threshold is None:
            threshold = cls.__DEFAULT_THRESHOLD

        np_err = np.geterr()
       
        np.seterr(divide='ignore', invalid='ignore')
    
        x = vars.shape[0]
        y = vars.shape[1]
   
        res_t = {}
        res_v = [{} for i in xrange(y)]
    
        if ui:
            ui.complete.value = 8 * (2 + num_features) # * (y + 1)
            ui.progress.value = 0
            ui.start()
    
        MI_t, H_t = Mrmr.__compute_mi(x, vars, targets, ui.progress if ui else None)
    
#         MI_v, H_v = np.zeros((y, y), dtype=float), np.zeros((y, y), dtype=float)
#     
#         for i in xrange(y):
#             MI_v[i, :], H_v[i, :] = compute_mi(x, vars, vars[:, i], ui.progress if ui else None)
    
        MIr_t = np.divide(MI_t, H_t)
    
        d_t = np.subtract(H_t, MI_t)
        D_t = np.divide(d_t, H_t)
    
#         MI_v = np.zeros((y, y), dtype=float) 
#         H_v = np.zeros((y, y), dtype=float)
#         d_v = np.zeros((y, y), dtype=float)
#         D_v = np.zeros((y, y), dtype=float)
#         for i in xrange(y):
#             for r in res_v[i].values():
#                 mi_v, h_v = r.get()
#                 MI_v[i, :] = np.add(MI_v[i, :], mi_v)
#                 H_v[i, :] = np.add(H_v[i, :], h_v)
#                 ui.progress.value += 1
#             d_v[i, :] = np.subtract(H_v[i, :], MI_v[i, :])
#             D_v[i, :] = np.divide(d_v[i, :], H_v[i, :])
    
        L_MI_t = MI_t.tolist()[0]
        L_MIr_t = MIr_t.tolist()[0]
    
#         L_d_t = d_t.tolist()[0]
#         L_D_t = D_t.tolist()[0]
    
        mi_vals = sorted(zip(range(len(L_MIr_t)), L_MIr_t), key=itemgetter(1), reverse=True) 
    
        idx, maxrel = mi_vals[0]
        mi_vars, h_vars, mir_vars = {}, {}, {}
       
        mi_vars[idx], h_vars[idx] = Mrmr.__compute_mi(x, vars, vars[:, idx], ui.progress if ui else None)
        mir_vars[idx] = np.divide(mi_vars[idx], h_vars[idx])
    
        # find related values 
        mu = sorted(zip(range(y), mir_vars[idx].tolist()[0]), key=itemgetter(1), reverse=True)
        related = [k for k in mu if k[1] > threshold and k[0] != idx]
        
        mrmr_vals = [(idx, maxrel, related)]
        mask_idxs = [idx]
    
        # do one extra because the sorting is sometimes off, do y-1 because we already include a feature by default 
        for k in xrange(min(num_features, y-1)):
            idx, maxrel, mrmr = sorted(
                [
                    (
                        idx,
                        maxrel,
                        # mRMR: MID then MIQ 
                        maxrel - np.sum(mir_vars[j][0, idx] for j, _, _ in mrmr_vals) / len(mrmr_vals) if method is MID else \
                        maxrel / np.sum(mir_vars[j][0, idx] for j, _, _ in mrmr_vals) / len(mrmr_vals)
                    ) \
                    for idx, maxrel in mi_vals[1:] if idx not in mask_idxs
                ], key=itemgetter(2), reverse=True)[0]
            mi_vars[idx], h_vars[idx] = Mrmr.__compute_mi(x, vars, vars[:, idx], ui.progress if ui else None)
            mir_vars[idx] = np.divide(mi_vars[idx], h_vars[idx])
    
            # find related values 
            mu = sorted(zip(range(y), mir_vars[idx].tolist()[0]), key=itemgetter(1), reverse=True)
            related = [k for k in mu if k[1] > cls.__DEFAULT_THRESHOLD and k[0] != idx]
    
            mrmr_vals.append((idx, mrmr, related))
            mask_idxs.append(idx)
        
        if ui:
            ui.join(0.1)
            if ui.is_alive():
                ui.terminate()
            stdout.write(' ' * 30 + '\r')
    
#         idx = mi_vals[0][0]
#         print 'I:', mi_vars[idx][0, idx]
#         print 'H:', h_vars[idx][0, idx]
#         print 'r:', mir_vars[idx][0, idx]
#         print 'd:', mi_vars[idx][0, idx] - h_vars[idx][0, idx]
#         print 'D:', (mi_vars[idx][0, idx] - h_vars[idx][0, idx]) / h_vars[idx][0, idx]
    
        # should be symmetric
#         assert(MI_v[0, 1] == MI_v[1, 0] and MI_v[0, y-1] == MI_v[y-1, 0])
#         assert(d_v[0, 1] == d_v[1, 0] and d_v[0, y-1] == d_v[y-1, 0])
#         assert(D_v[0, 1] == D_v[1, 0] and D_v[0, y-1] == D_v[y-1, 0])
        
#         print mi_t
#         print mi_v
    
        np.seterr(**np_err)
    
        return mi_vals[:num_features], sorted(mrmr_vals, key=itemgetter(1), reverse=True)[:num_features]
   
    def select(self, x, y):
        vars = np.matrix(x, dtype=bool)
        targets = np.matrix(y, dtype=bool)

        # transpose if necessary (likely if coming from array) 
        if targets.shape[0] == 1 and targets.shape[1] == vars.shape[0]:
            targets = targets.T
        elif targets.shape[1] != 1 or targets.shape[0] != vars.shape[0]:
            raise ValueError('`y\' should have as many entries as `x\' has rows.')

        self.__maxrel, self.__mrmr = None, None
        self.__maxrel, self.__mrmr = Mrmr.__mrmr_selection(self.num_features, self.method, targets, vars, self.threshold, ui)

        self.__colsize = vars.shape[1]
        self.__computed = True

    def features(self):
        if not self.__computed:
            raise StandardError('No mRMR model computed')

        if self.method is MAXREL:
            return [iv[0] for iv in self.__maxrel]
        else:
            return [ivr[0] for ivr in self.__mrmr]

    def subset(self, x):
        if not self.__computed:
            raise StandardError('No mRMR model computed')
        if x.shape[1] != self.__colsize:
            raise ValueError('model, number of features: shape mismatch')
        return x[:, Mrmr.features(self)]
