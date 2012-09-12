
from logging import getLogger
from math import copysign, sqrt
from operator import itemgetter
from re import compile as re_compile
from six import u
from unicodedata import combining

from mrmr import MRMR_LOGGER

from .._normalvalue import NormalValue


__all__ = ['IdepiResults']


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
        stat_prefixes[k] = sum([1 for c in u(k) if combining(c) == 0])

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


class IdepiResults(dict):

    def __init__(self, similar=True):
        self.__similar = similar

    def cv_results(self, cv_results, colnames):
        statsdict = cv_results.stats.todict()

        # remove minstat 'cause we don't want it here..
        if 'Minstat' in statsdict:
            del statsdict['Minstat']

        idxnames = {}
        weightsdict = {}

        for i in range(len(cv_results.extra)):
            # weights - 1 to account for bias term we added in LinearSvm.__bias
            assert(len(cv_results.extra[i]['features']) >= len(cv_results.extra[i]['weights']) - 1)
            # we only ever go up to the # of features selected by mRMR,
            # so the last weight (of the bias term) is ignored here, intentionally
            for j in range(len(cv_results.extra[i]['features'])):
                v = cv_results.extra[i]['weights'][j] if j < len(cv_results.extra[i]['weights']) else 0.
                k = cv_results.extra[i]['features'][j]
                r = set([idx for idx, _ in cv_results.extra[i]['similar'][k]]) if self.__similar else None
                name = colnames[k]
                idxnames[k] = name
                if name not in weightsdict:
                    weightsdict[name] = (NormalValue(int), set())
                weightsdict[name][0].append(int(copysign(1, v)))
                weightsdict[name][1].update(r)

        log = getLogger(MRMR_LOGGER)
        log.debug('mrmr index to name map: {%s}' % ', '.join(
            "%d: '%s'" % (
                idx, colnames[idx]
            ) for idx, name in sorted(idxnames.items(), key=itemgetter(0))
        ))

        numeric = re_compile(r'[^0-9]+')

        self['statistics'] = dict((k, { 'mean': v.mu, 'std': sqrt(v.sigma) }) for k, v in statsdict.items())
        self['weights'] = [{ 'position': k, 'value': { 'mean': v[0].mu, 'std': sqrt(v[0].sigma), 'N': len(v[0]) } } for k, v in sorted(
            weightsdict.items(),
            key=lambda x: int(numeric.sub('', x[0]))
        )]

        if self.__similar:
            for i, kv in enumerate(sorted(weightsdict.items(), key=lambda x: int(numeric.sub('', x[0])))):
                _, v = kv
                self['weights'][i]['similar'] = [colnames[j] for j in v[1]]

    def metadata(self, opts, N, balance, antibody, forward_select=None):
        self['metadata'] = {
            'sequences': N,
            'balance': balance,
            'features': opts.NUM_FEATURES if forward_select is None else forward_select,
            'discriminator': { 'orientation': 'gt', 'cutoff': opts.IC50 },
            'antibody': antibody,
            'folds': opts.CV_FOLDS
        }

    def predictions(self, ids, preds):
        assert len(ids) == len(preds), 'ids and preds are not the same length!'
        self['predictions'] = list(zip(ids, preds))

    def weights(self, weights):
        self['weights'] = weights

    def dumps(self, keys=None):
        if keys is None:
            keys = sorted(self.keys())
        ret = ['{\n']
        for i, key in enumerate(keys):
            if i > 0:
                ret.append(',\n')
            if key == 'metadata':
                ret.append(_dumps_metadata(self[key], 1))
            if key == 'statistics':
                ret.append(_dumps_statistics(self[key], 1))
            if key == 'weights':
                ret.append(_dumps_weights(self[key], 1, self.__similar))
            if key == 'predictions':
                ret.append(_dumps_predictions(self[key], 1))
        ret.append('\n}')
        return ''.join(ret)
