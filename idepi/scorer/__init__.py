
from warnings import catch_warnings, simplefilter

from numpy import (
    eye,
    mean,
    prod,
    seterr,
    zeros
    )

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

    def __call__(self, clf, X_pred, y_true):
        y_pred = clf.predict(X_pred)
        return Scorer.stats(y_true, y_pred)[self.__optstat]

    @staticmethod
    def __len__():
        return Scorer.NSTAT

    @staticmethod
    def __getitem__(key):
        return Scorer._NAMES[key]

    @staticmethod
    def stats(y_true, y_pred):
        vals = sorted(
            set(y_true.reshape((prod(y_true.shape),))) |
            set(y_pred.reshape((prod(y_pred.shape),)))
            )
        cm = confusion_matrix(y_true, y_pred, vals)
        # avoid division by zero errors
        def div(num, den):
            try:
                with catch_warnings(record=False):
                    simplefilter('ignore')
                    rv = num / den
            except ZeroDivisionError:
                rv = 0.0
            except FloatingPointError:
                rv = 0.0
            return rv
        def calc(tn, fp, fn, tp):
            npv = div(tn, tn + fn)
            ppv = div(tp, tp + fp) # precision
            sen = div(tp, tp + fn) # sensitivity / recall
            spe = div(tn, tn + fp) # specificity
            return npv, ppv, sen, spe
        nval = len(vals)
        eye_ = eye(nval)
        acc = div((eye_ * cm).sum(), cm.sum())
        if nval == 1:
            npv = ppv = sen = spe = 1.0
        elif nval == 2:
            npv, ppv, sen, spe = calc(*cm.reshape((4,)))
        else:
            npss = zeros((nval, 4))
            for i in range(nval):
                tp = cm[i, i]
                tn = (eye_ * cm).sum() - tp
                fp = cm[:, i].sum() - tp
                fn = cm[i, :].sum() - tp
                npss[i, :] = calc(tn, fp, fn, tp)
            npv, ppv, sen, spe = mean(npss, axis=0)
        f_1 = div(2 * ppv * sen, ppv + sen)
        stats = [None for _ in range(Scorer.NSTAT)]
        stats[Scorer.ACCURACY] = acc
        stats[Scorer.F1SCORE] = f_1
        stats[Scorer.MCC] = mcc(y_true, y_pred)
        stats[Scorer.NPV] = npv
        stats[Scorer.PPV] = ppv
        stats[Scorer.SENSITIVITY] = sen
        stats[Scorer.SPECIFICITY] = spe
        return stats
