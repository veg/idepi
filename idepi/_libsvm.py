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

from operator import itemgetter
from os import remove
from os.path import exists
from re import sub
from sys import stderr
from subprocess import Popen, PIPE
from tempfile import mkstemp

from numpy import mean, std

from _smldata import SmlData


class LibSVM(object):
    CSVC,     LINEAR      = 0
    NUSVC,    POLYNOMIAL  = 1
    ONECLASS, RBF         = 2
    EPSSVR,   SIGMOID     = 3
    NUSVR,    PRECOMPUTED = 4

    def __init__(self, train_data, test_data=None, model=None):
        # create a default trainfile is none exists
        self.model = model
        # set the other variables
        self.train_data = train_data
        self.test_data = test_data

    def __del__(self):
        if self._trainfile_is_temp:
            remove(self.trainfile)
        if self._testfile_is_temp:
            remove(self.testfile)
        if self._model_is_temp:
            remove(self.model)
        # pass

    @staticmethod
    def __parse_lsvm_output(opt_str, output):
        val = None
        stats = []
        for line in output.split("\n"):
            if line.find(opt_str) != -1:
                val = float(line.split()[-1][0:-1])
                stats = output.split("\n")[-9:-1]
        try:
            assert(val is not None and len(stats) > 0)
        except AssertionError, e:
            print >> stderr, output
            raise e
        return val, stats

    def grid_search(self, args=None, optimize="min2-5", c_begin=-5, c_end=15, c_step=2, nested=True, folds=5, svm=None):

        if not exists(self.trainfile) or self._trainfile_is_temp:
            if self.train_data is None or len(self.train_data) == 0:
                raise StandardError("ERROR: %s not found and no data exists to save there" % self.trainfile)
            if self.test_data is not None and len(self.test_data) == 0:
                raise StandardError("ERROR: test data is provided but empty!")
            LibSVM.save(self)
        
        if self.test_data is None and not nested:
            raise StandardError("ERROR: non-nested model-fitting requires test data")

        recip = 1
        if isinstance(c_step, float):
            recip = 1. / c_step
            c_begin, c_end = int(recip * c_begin), int(recip * c_end)
            c_step = 1
        c_range = [pow(2., 1. * c / recip) for c in xrange(c_begin, c_end + 1, c_step)]

        if optimize.lower() == "accuracy":
            opt_str = "Accuracy"
        elif optimize.lower() == "ppv":
            opt_str = "PPV"
        elif optimize.lower() == "npv":
            opt_str = "NPV"
        elif optimize.lower() == "sensitivity":
            opt_str = "Sensitivity"
        elif optimize.lower() == "specificity":
            opt_str = "Specificity"
        elif optimize.lower() == "f-score":
            opt_str = "F-score"
        elif optimize.lower() == "min2-5":
            opt_str = "min2-5"
        else:
            raise StandardError("ERROR: invalid optimization metric, valid metrics are:", \
              "accuracy, ppv, npv, sensitivity, specificity, f-score, and min2-5")

        svm_args = args if args is not None else []

        if not nested and "-v" in svm_args:
            raise StandardError("ERROR: You specified non-nested grid search but included a cross-validation parameter")
        elif nested:
            if "-v" in svm_args and folds is not None:
                raise StandardError("ERROR: you double-specified the number of cross-validation folds")
            elif "-v" not in svm_args and folds is None:
                raise StandardError("ERROR: you did not specify the number of cross-validation folds")
            else:
                svm_args += ["-v", "%d" % folds]

        best_c = 0.
        best_val = -1.
        best_stats = None
        best_weights = None

        for c in c_range:
            if nested:
                weights = None
                svm_process = Popen([self.svm_train] + svm_args + ["-c", "%s" % c, self.trainfile], close_fds = True, stdout = PIPE)
                output = svm_process.communicate()[0].strip()
            else:
                svm_train_process = Popen([self.svm_train] + svm_args + ["-c", "%s" % c, self.trainfile, self.model], close_fds = True, stdout = PIPE)
                svm_train_process.communicate()
                weights = LSVMModel(self.model, self.train_data.feature_names).weights()
                svm_predict_process = Popen([self.svm_predict] + [self.testfile, self.model, "/dev/null"], close_fds = True, stdout = PIPE)
                output = svm_predict_process.communicate()[0].strip()
            val, stats = LibSVM.__parse_lsvm_output(self, opt_str, output)
            if val > best_val:
                best_val = val
                best_c = c
                best_stats = stats 
                if not nested:
                    best_weights = weights

        if nested:
            # remove CV argument
            idx = svm_args.index("-v")
            svm_args.pop(idx)
            svm_args.pop(idx)
            # train the damn thing
            svm_train_process = Popen([self.svm_train] + svm_args + ["-c", "%s" % best_c, self.trainfile, self.model], close_fds = True, stdout = PIPE)
            svm_train_process.communicate()
            best_weights = LSVMModel(self.model, self.train_data.feature_names).weights()
            if self.test_data is not None:
                svm_predict_process = Popen([self.svm_predict] + [self.testfile, self.model, "/dev/null"], close_fds = True, stdout = PIPE)
                output = svm_predict_process.communicate()[0].strip()
                _, best_stats = LibSVM.__parse_lsvm_output(self, opt_str, output)

        if best_c == 0.:
            print >> stderr, "Unable to find best regularization parameter!"
            exit(-1)

        stats_dict = {}
        for line in stats:
            if '=' not in line:
                continue
            k, v = line.split('=')
            k = sub(r"^\d\)\s+Cross\s+Validation\s+", "", k)
            stats_dict[k] = float(v.rstrip('%'))

        return best_c, stats_dict, best_weights

    def save(self):
        fh = open(self.trainfile, 'w')
        for r in self.train_data:
            # LIBSVM format "assumes" model features are 1-indexed, so use (f[0] + 1, f[1]) below
            print >> fh, "%d" % r.value, ' '.join(["%d:%s" % (f[0] + 1, f[1]) for f in sorted(r.features.items(), key = itemgetter(0)) if f[1] != 0])
        fh.close()
        if self.test_data is not None:
            fh = open(self.testfile, 'w')
            for r in self.test_data:
                print >> fh, "%d" % r.value, ' '.join(["%d:%s" % (f[0] + 1, f[1]) for f in sorted(r.features.items(), key = itemgetter(0)) if f[1] != 0])
            fh.close()
