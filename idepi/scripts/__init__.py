
from ._discrete import main as discrete
from ._learn import main as learn
from ._phylo import main as phylo
from ._predict import main as predict
# from ._regressor import main as regressor
from ._sto2fa import main as sto2fa
from ._tree import main as tree

__all__ = [
    'discrete',
    'learn',
    'phylo',
    'predict',
#     'regressor',
    'sto2fa',
    'tree'
]
