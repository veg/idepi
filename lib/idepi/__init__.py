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

__version__ = '0.6'

from . import _phylogzfile as PhyloGzFile

from ._alphabet import *
from ._classextractor import *
from ._common import *
from ._discretizer import *
from ._hmmer import *
from ._hyphy import *
from ._linearsvm import *
from ._logging import *
from ._mlpy import *
from ._naivefilter import *
from ._normalvalue import *
from ._orflist import *
from ._phylo import *
from ._phylofilter import *
from ._randomsequences import *
# from _simulatedepitope import *
from ._simulation import *
from ._smldata import *
from ._sparsepartitioning import *
from ._util import *


__all__ = ['IDEPI_LOGGER', 'PhyloGzFile']
__all__ += _alphabet.__all__
__all__ += _classextractor.__all__
__all__ += _common.__all__
__all__ += _discretizer.__all__
__all__ += _hmmer.__all__
__all__ += _hyphy.__all__
__all__ += _linearsvm.__all__
# don't import logging here, we just want to setup
__all__ += _mlpy.__all__
__all__ += _naivefilter.__all__
__all__ += _normalvalue.__all__
__all__ += _orflist.__all__
__all__ += _phylo.__all__
__all__ += _phylofilter.__all__
__all__ += _randomsequences.__all__
# __all__ += _simulatedepitope.__all__
__all__ += _simulation.__all__
__all__ += _smldata.__all__
__all__ += _sparsepartitioning.__all__
__all__ += _util.__all__

try:
    from ._logo import *
    __all__ += _logo.__all__
except ImportError:
    pass

_setup_log()
