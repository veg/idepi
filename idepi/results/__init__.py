
from operator import itemgetter
from re import compile as re_compile
from six import u
from unicodedata import combining

import numpy as np

from ..normalvalue import NormalValue


__all__ = ['Results']


def _dumps_metadata(meta, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix + '"metadata": {\n'

    name_len = max(len(k) for k in meta.keys()) + 3
    output = (prefix + '  %-*s %s' % (
        name_len,
        '"%s":' % k,
        '"%s"' % v if isinstance(v, str) else \
        ' { %s }' % ', '.join(
            ('"%s": %s' % (
                k,
                '"%s"' % v if isinstance(v, str) else
                '%.6g' % v if isinstance(v, float) else
                '%s' % str(v)
            ) for k, v in v.items())
        ) if isinstance(v, dict) else \
        ' %.6g' % v if isinstance(v, float) else \
        ' %s' % str(v)
    ) for k, v in sorted(meta.items(), key=itemgetter(0)))

    return buf + ',\n'.join(output) + '\n' + prefix + '}'


def _dumps_predictions(preds, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix + '"predictions": {\n'

    id_len = max(len(id) for id, _ in preds) + 3

    output = (prefix + '  %-*s %d' % (
        id_len,
        '"%s":' % id,
        pred
    ) for id, pred in preds)

    return buf + ',\n'.join(output) + '\n' + prefix + '}'


def _dumps_statistics(stats, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix

    buf += '"statistics": {\n'

    stat_prefixes = {}
    for k in stats.keys():
        stat_prefixes[k] = sum(1 for c in u(k) if combining(c) == 0)

    stat_len = max(stat_prefixes.values())
    mean_len = max(len('%.6f' % v['mean']) for v in stats.values())
    std_len = max(len('%.6f' % v['std']) for v in stats.values())
    fmt = '{ "mean": %%%d.6f, "std": %%%d.6f }' % (mean_len, std_len)
    output = (prefix + '  %s%s %s' % (
        '"%s":' % k,
        ' ' * (stat_len - stat_prefixes[k]),
        fmt % (v['mean'], v['std'])
    ) for k, v in sorted(stats.items(), key=itemgetter(0)))

    return buf + ',\n'.join(output) + '\n' + prefix + '}'


def _dumps_weights(weights, ident=0, similar=True):
    numeric = re_compile(r'[^0-9]+')
    def weightkey(v):
        return int(numeric.sub('', v['position']))

    prefix = ' ' * 2 * ident

    buf = prefix + '"weights": [\n'

    if len(weights) > 0:
        similar = False if 'similar' not in weights[0] else similar
        name_len = max(len(v['position']) for v in weights) + 3
        if isinstance(weights[0]['value'], dict):
            mean_len = max(len('% .6f' % v['value']['mean']) for v in weights)
            std_len = max(len('%.6f' % v['value']['std']) for v in weights)
            N_len = max(len('%d' % v['value']['N']) for v in weights)
            fmt = '{ "mean": %%%d.6f, "std": %%%d.6f, "N": %%%dd }' % (mean_len, std_len, N_len)
        elif isinstance(weights[0]['value'], int):
            val_len = max(len('% d' % v['value']) for v in weights)
            fmt = '%% %dd' % val_len
        else:
            raise RuntimeError('someone is fucking with us')
        if similar:
            similar_len = max(len(', '.join('"%s"' % r for r in v['similar'])) for v in weights)
            output = (prefix + '  { "position": %-*s "value": %s, "similar": [ %-*s ] }' % (
                name_len, '"%s",' % v['position'],
                fmt % (
                    (
                        v['value']['mean'],
                        v['value']['std'],
                        v['value']['N']
                    ) if isinstance(v['value'], dict) else (
                        v['value']
                    )
                ),
                similar_len, ', '.join('"%s"' % r for r in sorted(v['similar'], key=lambda v: int(numeric.sub('', v))))
            ) for v in sorted(weights, key=weightkey))
        else:
            output = (prefix + '  { "position": %-*s "value": %s }' % (
                name_len, '"%s",' % v['position'],
                fmt % (
                    (
                        v['value']['mean'],
                        v['value']['std'],
                        v['value']['N']
                    ) if isinstance(v['value'], dict) else (
                        v['value']
                    )
                )
            ) for v in sorted(weights, key=weightkey))

    return buf + ',\n'.join(output) + '\n' + prefix + ']'


class Results(dict):

    def __init__(self, labels, scorer, similar=False):
        super(Results, self).__init__()
        self.__labels = labels
        self.__similar = similar
        self.__nfeat = NormalValue(int)
        self.__nfold = 0
        self.__npos = 0
        self.__ntotal = 0
        self.__scorer = scorer
        self.__stats = [
            NormalValue(float, name=self.__scorer[i])
            for i in range(len(self.__scorer))
            ]
        self.__coefs = [NormalValue(int) for _ in range(len(labels))]
        self.__valid = False

    def add(self, y_true, y_pred, coefs):
        self.__nfeat.add(len(coefs))
        self.__nfold += 1
        self.__npos += (y_true > 0).sum()
        self.__ntotal += np.prod(y_true.shape)
        for i, coef in coefs.items():
            self.__coefs[i].add(coef)
        # update stats
        stats_ = self.__scorer.stats(y_true, y_pred)
        for i in range(len(self.__scorer)):
            self.__stats[i].add(stats_[i])
        # invalidate state
        self.__valid = False

    @property
    def optstat(self):
        return self.__stats[self.__scorer.optstat]

    def __lt__(self, other):
        return self.optstat < other.optstat

    def __le__(self, other):
        return self.optstat <= other.optstat

    def __eq__(self, other):
        return self.optstat == other.optstat

    def __ne__(self, other):
        return self.opstat != other.optstat

    def __gt__(self, other):
        return self.optstat > other.optstat

    def __ge__(self, other):
        return self.optstat >= other.optstat

    def __getitem__(self, key):
        Results.__compute(self)
        return super(Results, self).__getitem__(key)

    def __compute(self):
        if self.__valid:
            return

        self['metadata'] = {
            'sequences': self.__ntotal,
            'balance': self.__npos / self.__ntotal,
            'features': int(self.__nfeat.mean),
            'folds': self.__nfold
        }

        self['statistics'] = dict(
            (
                stat.name,
                dict(mean=stat.mean, std=stat.std)
                )
            for stat in self.__stats
            )

        weights = []
        for i, label in enumerate(self.__labels):
            v = self.__coefs[i]
            if len(v) == 0:
                continue
            weights.append(
                dict(
                    position=label,
                    value=dict(
                        mean=v.mean,
                        std=v.std,
                        N=len(v)
                        )
                    )
                )
        self['weights'] = weights
        self.__valid = True

    def metadata(self, antibody, ic50):
        Results.__compute(self)
        self['metadata']['antibody'] = antibody
        self['metadata']['discriminator'] = { 'orientation': 'gt', 'cutoff': ic50 }

    def predictions(self, ids, preds):
        assert len(ids) == len(preds), 'ids and preds are not the same length!'
        self['predictions'] = list(zip(ids, preds))

    def dumps(self, keys=None):
        Results.__compute(self)
        if keys is None:
            keys = sorted(self.keys())
        ret = ['{\n']
        for i, key in enumerate(keys):
            if i > 0:
                ret.append(',\n')
            if key == 'metadata':
                ret.append(_dumps_metadata(self[key], 1))
            elif key == 'statistics':
                ret.append(_dumps_statistics(self[key], 1))
            elif key == 'weights':
                ret.append(_dumps_weights(self[key], 1, self.__similar))
            elif key == 'predictions':
                ret.append(_dumps_predictions(self[key], 1))
            else:
                raise RuntimeError("unknown key: '{0}'".format(key))
        ret.append('\n}')
        return ''.join(ret)
