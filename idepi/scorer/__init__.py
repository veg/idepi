
from numpy import seterr

from sklearn.metrics import (
    confusion_matrix,
    matthews_corrcoef
    )


__all__ = ['Scorer', 'mcc']


def mcc(*args, **kwargs):
    settings = seterr(invalid='ignore')
    mcc_ = matthews_corrcoef(*args, **kwargs)
    seterr(**settings)
    return mcc_

class Scorer:
    ACCURACY, F1SCORE, MCC, NPV, PPV, SENSITIVITY, SPECIFICITY, NSTAT = range(8)
    _NAMES = ['accuracy', 'f1score', 'mcc', 'npv', 'ppv', 'sensitivity', 'specificity']

    def __init__(self, optstat=None):
        if optstat is None:
            optstat = Scorer.MCC
        else:
            if optstat not in range(Scorer.NSTAT):
                raise ValueError(
                    'optstat must be one of Scorer.{ACCURACY, F1SCORE, MCC, NPV, PPV, SENSITIVITY, SPECIFICITY}'
                    )
        self.__optstat = optstat

    @property
    def optstat(self):
        return self.__optstat

    def __call__(self, y_true, y_pred):
        return Scorer.stats(y_true, y_pred)[self.__optstat]

    @staticmethod
    def __len__():
        return Scorer.NSTAT

    @staticmethod
    def __getitem__(key):
        return Scorer._NAMES[key]

    @staticmethod
    def stats(y_true, y_pred):
        tp, fn, fp, tn = confusion_matrix(y_true, y_pred, [True, False]).reshape((4,))
        acc = (tp + tn) / (tp + tn + fp + fn)
        npv = tn / (tn + fn)
        ppv = tp / (tp + fp) # precision
        sen = tp / (tp + fn) # sensitivity / recall
        spe = tn / (tn + fp) # specificity
        f_1 = 2 * ppv * sen / (ppv + sen)
        stats = [None for _ in range(Scorer.NSTAT)]
        stats[Scorer.ACCURACY] = acc
        stats[Scorer.F1SCORE] = f_1
        stats[Scorer.MCC] = mcc(y_true, y_pred)
        stats[Scorer.NPV] = npv
        stats[Scorer.PPV] = ppv
        stats[Scorer.SENSITIVITY] = sen
        stats[Scorer.SPECIFICITY] = spe
        return stats
