
from copy import copy
from re import _pattern_type as RegexObject


__all__ = [
    'RegexGlobber'
    ]


class RegexGlobber(dict):

    def __init__(self, regexes, alignment=None, refidx=None):
        if not all(isinstance(regex, RegexObject) for regex in regexes):
            raise ValueError(
                "regexes must be an iterable of compiled RegexObjects"
                )
        self._regexes = list(regexes)
        self._rows = {}
        if alignment is not None:
            self(alignment, refidx)

    def reset(self):
        self._rows.clear()

    def __call__(self, alignment, refidx=None):
        r = 0
        for i, seq in enumerate(alignment):
            if i == refidx:
                continue
            m = None
            for regex in self._regexes:
                m = regex.match(seq.id)
                if m is not None:
                    break
            if m is None:
                raise ValueError(
                    "key '%s' does not match supplied regular expression" % seq.id
                    )
            else:
                row, _ = RegexGlobber.__row_weight(m)
                if row not in self._rows:
                    self._rows[row] = r
                    r += 1

    @staticmethod
    def __row_weight(m):
        groups = m.groups()
        try:
            weight = m.group('weight')[0]
            groups = [g for g in groups if g != weight]
            weight = float(weight)
        except IndexError:
            weight = 1
        row = ''.join(groups)
        return row, weight

    def __getitem__(self, seqid):
        m = None
        for regex in self._regexes:
            m = regex.match(seqid)
            if m is not None:
                break
        if m is None:
            raise ValueError(
                "key '%s' does not match supplied regular expression" % seqid
                )
        else:
            row, weight = RegexGlobber.__row_weight(m)
            return self._rows[row], weight

    def __len__(self):
        return len(self._rows)

    def __setitem__(self, key, value):
        if key in self._rows:
            raise ValueError('Globbers cannot be updated')
        else:
            self._rows[key] = value

    def __add__(self, other):
        new = RegexGlobber(self._regexes)
        new._rows = copy(self._rows)
        new += other
        return new

    def __iadd__(self, other):
        self._regexes.extend(other._regexes)
        for row, idx in other._rows.items():
            self[row] = len(self) + idx
        return self
