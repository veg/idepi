
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
    output = (
        prefix + '  %-*s %s' % (
            name_len,
            '"%s":' % k,
            '"%s"' % v if isinstance(v, str) else
            ' { %s }' % ', '.join((
                '"%s": %s' % (
                    k_,
                    '"%s"' % v_ if isinstance(v_, str) else
                    '%.6g' % v_ if isinstance(v_, float) else
                    '%s' % str(v_)
                    ) for k_, v_ in v.items())
            ) if isinstance(v, dict) else
            ' [ %s ]' % ', '.join((
                '"%s"' % v_ if isinstance(v_, str) else
                '%.6g' % v_ if isinstance(v_, float) else
                '%s' % str(v_)
                ) for v_ in v
            ) if isinstance(v, (list, tuple)) else
            ' %.6g' % v if isinstance(v, float) else
            ' %s' % str(v)
        )
        for k, v in sorted(meta.items(), key=itemgetter(0))
        if not isinstance(v, int) or v
        )

    return buf + ',\n'.join(output) + '\n' + prefix + '}'


def _dumps_predictions(preds, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix + '"predictions": {\n'

    id_len = max(len(id) for id, _ in preds) + 3

    output = (prefix + '  %-*s % d' % (
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

    if len(weights) == 0:
        output = ''
    else:
        similar = False if 'similar' not in weights[0] else similar
        name_len = max(len(v['position']) for v in weights) + 3
        N_len = max(len('%d' % v['N']) for v in weights) + 1
        # rank
        mean_len = max(len('%.2f' % v['rank']['mean']) for v in weights)
        std_len = max(len('%.2f' % v['rank']['std']) for v in weights)
        rank_fmt = '{ "mean": %%%d.2f, "std": %%%d.2f }' % (mean_len, std_len)
        if isinstance(weights[0]['value'], dict):
            mean_len = max(len('% .2f' % v['value']['mean']) for v in weights)
            std_len = max(len('%.2f' % v['value']['std']) for v in weights)
            val_fmt = '{ "mean": %% %d.2f, "std": %%%d.2f }' % (mean_len, std_len)
        elif isinstance(weights[0]['value'], int):
            val_len = max(len('% d' % v['value']) for v in weights)
            val_fmt = '%% %dd' % val_len
        else:
            raise RuntimeError('someone is fucking with us')
        if similar:
            similar_len = max(len(', '.join('"%s"' % r for r in v['similar'])) for v in weights)
            output = ',\n'.join(prefix + '  { "position": %-*s "N": %-*s "rank": %s, "value": %s, "similar": [ %-*s ] }' % (
                name_len, '"%s",' % v['position'],
                N_len, '%d,' % v['N'],
                rank_fmt % (
                    v['rank']['mean'],
                    v['rank']['std']
                    ),
                val_fmt % (
                    (
                        v['value']['mean'],
                        v['value']['std']
                        ) if isinstance(v['value'], dict) else (
                            v['value']
                        )
                    ),
                similar_len,
                ', '.join(
                    '"%s"' % r for r in sorted(
                        v['similar'],
                        key=lambda v: int(numeric.sub('', v))
                        )
                    )
                ) for v in sorted(weights, key=weightkey)) + '\n'
        else:
            output = ',\n'.join(prefix + '  { "position": %-*s "N": %-*s "rank": %s, "value": %s }' % (
                name_len, '"%s",' % v['position'],
                N_len, '%d,' % v['N'],
                rank_fmt % (
                    v['rank']['mean'],
                    v['rank']['std']
                    ),
                val_fmt % (
                    (
                        v['value']['mean'],
                        v['value']['std']
                        ) if isinstance(v['value'], dict) else (
                            v['value']
                        )
                    )
                ) for v in sorted(weights, key=weightkey)) + '\n'

    return buf + output + prefix + ']'


class Results(dict):

    def __init__(self, labels, scorer, similar=0.0):
        super(Results, self).__init__()
        self.__labels = labels
        self.__similar = similar > 0.0
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
        self.__ranks = [NormalValue(int) for _ in range(len(labels))]
        self.__valid = False

    def add(self, y_true, y_pred, coefs, ranks):
        # skip fold data if not present
        if y_true is not None and y_pred is not None:
            self.__nfold += 1
            stats_ = self.__scorer.stats(y_true, y_pred)
            for i in range(len(self.__scorer)):
                self.__stats[i].add(stats_[i])

        if y_true is not None:
            self.__npos += (y_true > 0).sum()
            self.__ntotal += np.prod(y_true.shape)

        self.__nfeat.add(len(coefs))

        for i, coef in coefs.items():
            self.__coefs[i].add(coef)
        for i, rank in ranks.items():
            self.__ranks[i].add(rank)

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

        # avoid division by zero errors
        balance = (
            self.__npos / self.__ntotal
            if self.__ntotal else 0.0
            )

        self['metadata'] = {
            'sequences': self.__ntotal,
            'balance': balance,
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
            r = self.__ranks[i]
            if len(v) == 0:
                continue
            weights.append(
                dict(
                    position=label,
                    N=len(v),
                    rank=dict(
                        mean=r.mean,
                        std=r.std
                        ),
                    value=dict(
                        mean=v.mean,
                        std=v.std,
                        )
                    )
                )
        self['weights'] = weights
        self.__valid = True

    def metadata(self, antibodies, label):
        Results.__compute(self)
        self['metadata']['antibodies'] = antibodies
        self['metadata']['label'] = label

    def predictions(self, ids, preds):
        assert len(ids) == len(preds), 'ids and preds are not the same length!'
        self['predictions'] = list(zip(ids, preds))

    def dumps(self, keys=None):
        Results.__compute(self)
        if keys is None:
            keys = sorted(self.keys() if self.__nfold else set(self.keys()) - set(['statistics']))
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
