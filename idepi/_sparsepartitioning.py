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
from os import chdir, getcwd, listdir, remove, rmdir
from os.path import exists, isfile, join
from subprocess import PIPE, Popen
from sys import exit, stderr
from tempfile import mkdtemp
from types import IntType, FloatType

from _smldata import SmlData


__all__ = ['SpDeterministicFeature', 'SpDeterministicModel', 'SpMcmcFeature', 'SpMcmcPair', 'SpMcmcModel', 'SparsePartitioning'] 


class SpDeterministicFeature(object):

    def __init__(self, idx, name):
        self.idx = idx
        self.name = name

    def __eq__(self, other):
        if self.idx == other.idx and self.name == other.name:
            return True
        return False


class SpDeterministicModel(object):

    def __init__(self, response, features, assocnum=None, maxvar=None):
        self.assocnum = assocnum
        self.response = response
        self.features = features
        self.maxvar = maxvar

    def __eq__(self, other):
        try:
            assert(len(self.features) == len(other.features))
            A = sorted(self.features, key=lambda x: x.idx)
            B = sorted(other.features, key=lambda x: x.idx)
            for i in xrange(len(A)):
                assert(A[i] == B[i])
            return True
        except AssertionError:
            return False


class SpMcmcFeature(object):

    def __init__(self, idx, name, post):
        self.idx = idx
        self.name = name
        self.post = post


class SpMcmcPair(object):

    def __init__(self, features, post, response):
        self.features = features
        self.post = post
        self.response = response


class SpMcmcModel(object):

    def __init__(self, features, pairs):
        self.features = features
        self.pairs = pairs


class SparsePartitioning(object):
    __DEFAULT_SP_ITER = 250
    __DEFAULT_SP_DET = 'spdet'
    __DEFAULT_SP_MCMC = 'sp'
    __DEFAULT_SP_METHOD = 'det'
    # 1: linear model (continuous response)
    # 2: logistic model (0/1 responses)
    # 3: probit model (0/1 responses), uses 'Albert Chib Latent Variable Representation'
    __DEFAULT_SP_OPTIONS = dict(iterno=__DEFAULT_SP_ITER, modtype=2)

    __KNOWN_SP_FILES = (
        'bp.txt',
        'chr.txt',
        'comadj.txt',
        'combells.txt',
        'comcor.txt',
        'comkeep.txt',
        'comrev.txt',
        'comtrue.txt',
        'data.txt',
        'detpost.txt',
        'detres.txt',
        'detsteps.txt',
        'dettop.txt',
        'ematrix.txt',
        'input.txt',
        'mod.txt',
        'phen.txt',
        'prior.txt',
        'spanova.txt',
        'spbf.txt',
        'spcore.txt',
        'spcoreb.txt',
        'spcv.txt',
        'spinter.txt',
        'spinterall.txt',
        'spiter.txt',
        'sppart.txt',
        'sppartb.txt',
        'sppost.txt',
        'spres.txt',
        'spshared.txt',
        'spsize.txt'
    )
    
    __POSITIVE_INTEGER = lambda x: type(x) is IntType and x > 0
    __POSITIVE_FLOAT = lambda x: type(x) in (IntType, FloatType) and x > 0.

    __KNOWN_SP_OPTIONS = dict(
        oldn=(True, __POSITIVE_INTEGER),
        oldN=(True, __POSITIVE_INTEGER),
        iterno=(True, __POSITIVE_INTEGER),
        modtype=(True, lambda x: x in (1, 2, 3)),
        M=(False, None),
        modenter=(False, None),
        priorent=(False, None),
        pprior=(None, lambda x: type(x) in (IntType, FloatType)),
        priorass=(None, lambda x: type(x) in (IntType, FloatType)),
        maxT=(None, __POSITIVE_INTEGER),
        maxD=(None, __POSITIVE_INTEGER),
        copies=(None, __POSITIVE_INTEGER),
        missingvalue=(False, None),
        regvar=(None, __POSITIVE_FLOAT),
        split=(None, lambda x: x in (0, 1)),
        sskip=(None, __POSITIVE_INTEGER),
        supdate=(None, __POSITIVE_INTEGER),
        tol=(None, lambda x: type(x) in (IntType, FloatType)),
        randomise=(None, lambda x: x in (0, 1)),
        maxiter=(None, __POSITIVE_INTEGER),
        numcor=(None, __POSITIVE_INTEGER),
        corvar=(None, __POSITIVE_FLOAT),
        posenter=(False, None),
        ldwind=(None, __POSITIVE_INTEGER),
        ldmax=(None, __POSITIVE_FLOAT)
    )

    def __init__(self, data, feature_names, dirname=None):
        if dirname is None:
            dirname = mkdtemp()
            self.__dir_is_temp = True
        else:
            self.__dir_is_temp = False
        self.dirname = dirname
        # must be data.txt, input.txt, phen.txt as required by SP
        self.filenames = dict(data=join(self.dirname, 'data.txt'),
                              input=join(self.dirname, 'input.txt'),
                              phen=join(self.dirname, 'phen.txt'))
        self.feature_names = feature_names
        self.data = data
        try:
            assert(self.data.feature_names == self.feature_names)
        except AssertionError, e:
            print >> stderr, 'ERROR: feature names do not match those of the SmlData, aborting'
            raise e
        self.method = None
        self.__options = None
        self.output = None

    def __del__(self):
        if self.__dir_is_temp:
            del_ = True
            for file_ in listdir(self.dirname):
                filepath = join(self.dirname, file_)
                if isfile(filepath):
                    if file in self.__KNOWN_SP_FILES:
                        remove(filepath)
                    else:
                        print >> stderr, 'WARNING: will not remove temporary directory %s, for it contains an unknown file %s' % (self.dirname, file_)
                        del_ = False
                else:
                    print >> stderr, 'WARNING: will not remove temporary directory %s, for it contains an unknown directory %s' % (self.dirname, file_)
                    del_ = False
            if del_:
                rmdir(self.dirname)

    # do this after error checking or you'll unintentionally set your options upon failure
    def __merge_options(self, options):
        for k, v in options.items():
            if k not in self.__KNOWN_SP_OPTIONS:
                print >> stderr, 'WARNING: ignoring option %s, for it isn\'t a valid option' % k
                del options[k]
            elif self.__KNOWN_SP_OPTIONS[k][0] is False:
                print >> stderr, 'WARNING: ignoring option %s, for it cannot be set by this module' % k
            elif self.__KNOWN_SP_OPTIONS[k][1](v) is False:
                print >> stderr, 'WARNING: ignoring option %s, for %s isn\'t valid for this option' % (k, str(v))
            else:
                self.__options[k] = v

    def __validate_options(self):
        # iterate over all known options, and if it's required, make sure it's both there and valid
        for k, v in self.__KNOWN_SP_OPTIONS.items():
            if v[0]: # the options is True, as in not in (False, None)
                if k not in self.__options:
                    print >> stderr, 'ERROR: option %s must be set!' % k
                    raise StandardError
                elif v[1](self.__options[k]) is False:
                    print >> stderr, 'ERROR: option %s has invalid value %s' % (k, str(self.__options[k]))
                    raise StandardError
        # iterate over all set options, and make sure it's valid
        for k, v in self.__options.items():
            if self.__KNOWN_SP_OPTIONS[k][1](v) is False:
                print >> stderr, 'ERROR: option %s has invalid value %s' % (k, v)
                raise StandardError

    def run(self, iter=None, method=None, sp=None, options=None, fixdet=True):
        if iter is None:
            iter = self.__DEFAULT_SP_ITER
        if method is None:
            method = self.__DEFAULT_SP_METHOD
        if options is None:
            options = self.__DEFAULT_SP_OPTIONS
        if method not in ('det', 'mcmc'):
            print >> stderr, 'ERROR: method must be either deterministic (\'det\') or monte carlo markov chain (\'mcmc\')'
            raise StandardError
        if sp is None:
            sp = self.__DEFAULT_SP_DET if method == 'det' else self.__DEFAULT_SP_MCMC
        # this is a convenience bit for appropriately transforming the function into its deterministic cousin
        # if we want to use the deterministic method, and is controlled by the fixdet parameter
        if fixdet and method == 'det':
            if match(r'sp$', sp):
                sp += 'det'
        save = False
        for reqfile in self.filenames.values():
            if not exists(reqfile) or self.__dir_is_temp:
                if self.feature_names is None or len(self.feature_names) == 0 or \
                    self.data is None or len(self.data) == 0:
                    print >> stderr, 'ERROR: %s not found and/or no data exists to save there' % reqfile
                    raise StandardError
                save = True
        if save:
            SparsePartitioning.save(self, options)
        cwd = getcwd()
        chdir(self.dirname)
        proc = Popen([sp, iter], close_fds=True, stdout=PIPE)
        self.output = proc.communicate()[0]
        chdir(cwd)
        SparsePartitioning.parse(self)
        self.method = method
        self.model

    def save(self, options=None):
        if options is None:
            options = self.__DEFAULT_SP_OPTIONS
        # set, validate, and print some configuration options before we write files
        oldn, oldN = len(self.data), len(self.feature_names)
        SparsePartitioning.__merge_options(self, dict(oldn=oldn, oldN=oldN))
        floats = False
        for r in self.data:
            if type(r.value) is FloatType:
                floats = True
        if self.__options is None:
            SparsePartitioning.__merge_options(self, options)
        if self.__options['modtype'] in (2, 3) and floats:
            print >> stderr, 'ERROR: modtype %i requires binary, not continuous, responses' % self.__options['modtype']
            raise StandardError
        SparsePartitioning.__validate_options(self)
        fh = open(self.filenames['input'], 'w')
        for option, value in self.__options.items():
            if value is not None and value is not False:
                print >> fh, '%s=%s' % (option, str(value))
        fh.close()
        # have to convert to column format
        preds = dict([(i, {}) for i in xrange(oldN)])
        for i in xrange(len(self.data)):
            for p, v in self.data[i].features.items():
                preds[p][i] = v
        # print the predictors (genotypes?)
        fh = open(self.filenames['data'], 'w')
        for p in sorted(preds.keys()):
            print >> fh, ' '.join([1 if preds[p][i] else 0 for i in xrange(oldn)])
        fh.close()
        # print the phenotypes
        fh = open(self.filenames['phen'], 'w')
        for r in self.data:
            print >> fh, r.value
        fh.close()

    def parse(self):
        finalmodel = None
        if self.method == 'det':
            models = []
            fh = open(join(self.dirname, 'detres.txt'))
            for line in fh:
                line = line.strip()
                response, predstr = line.split(' ', 1)
                response = int(response)
                predidxs = predstr.split(' ')
                # use a -1 to 0-index things
                predidxs = [int(p) - 1 for p in predidxs]
                features = [SpDeterministicFeature(i, self.feature_names[i]) for i in predidxs]
                lastmodel = SpDeterministicModel(response, features) 
            fh.close()
            fh = open(join(self.dirname, 'detsteps.txt'))
            for line in fh:
                line = line.strip()
                response, assocnum, maxvar, predstr = line.split(' ', 3)
                response, assocnum, maxvar = int(response), int(assocnum), float(maxvar)
                predidxs = predstr.split(' ')
                predidxs = [int(p) - 1 for p in predidxs]
                features = [SpDeterministicFeature(i, self.feature_names[i]) for i in predidxs]
                models.append(SpDeterministicModel(response, features, assocnum, maxvar))
            fh.close()
            for model in models:
                if model == finalmodel:
                    finalmodel = model
        else:
            fh = open(join(self.dirname, 'spres.txt'))
            features = []
            c = 0
            for line in fh:
                line = line.strip()
                post, _ = line.split(' ')
                features.append(SpMcmcFeature(c, self.feature_names[c], float(post)))
                c += 1
            fh.close()
            fh = open(join(self.dirname, 'spinter.txt'))
            pairs = []
            for line in fh:
                line = line.strip()
                idx1, idx2, post, response = line.split(' ')
                idx1, idx2, post, response = int(idx1), int(idx2), float(post), int(response)
                pairs.append(SpMcmcPair([features[idx1-1], features[idx2-1]], post, response))
            fh.close()
            finalmodel = SpMcmcModel(sorted(features, key=lambda x: x.post, reverse=True),
                                     sorted(pairs, key=lambda x: x.post, reverse=True)) 
        self.model = finalmodel
        return self.model
