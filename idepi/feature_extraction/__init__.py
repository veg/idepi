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

from sklearn.pipeline import FeatureUnion as FeatureUnion_

from idepi.feature_extraction._msavectorizer import *
from idepi.feature_extraction._msavectorizerpairwise import *
from idepi.feature_extraction._msavectorizerregex import *
from idepi.feature_extraction._msavectorizerregexpairwise import *

__all__ = ['FeatureUnion']
__all__ += _msavectorizer.__all__
__all__ += _msavectorizerpairwise.__all__
__all__ += _msavectorizerregex.__all__
__all__ += _msavectorizerregexpairwise.__all__


class FeatureUnion(FeatureUnion_):
    def get_feature_names(self):
        feature_names = []
        for _, trans in self.transformer_list:
            feature_names.extend(trans.get_feature_names())
        return feature_names
