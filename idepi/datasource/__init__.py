
from __future__ import division, print_function

from csv import reader as csv_reader, Sniffer as csv_sniffer
from json import dumps as json_dumps
from os.path import basename, splitext
from re import compile as re_compile
from sqlite3 import connect
from textwrap import dedent
from warnings import warn

from Bio import SeqIO
from Bio.Alphabet import Gapped, generic_nucleotide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt.orflist import OrfList


__all__ = ['DataSource']


def DataSource(s):
    args = s.split(',')
    if len(args) == 1:
        return Sqlite3Db(*args)
    if len(args) == 2:
        return MonogramData(*args)
    else:
        raise ValueError("invalid data specification: '%s'" % s)


class Sqlite3Db:

    def __init__(self, filename):
        self.__filename = filename

    @property
    def basename_root(self):
        return splitext(basename(self.__filename))[0]

    @property
    def filename(self):
        return self.__filename

    @property
    def antibodies(self):
        conn = connect(self.__filename)
        cur = conn.cursor()
        cur.execute('select distinct ANTIBODY from ANTIBODY')
        valid_antibodies = [r[0] for r in cur]
        conn.close()
        return valid_antibodies

    @property
    def labels(self):
        conn = connect(self.__filename)
        cur = conn.cursor()
        cur.execute('select distinct TYPE from NEUT_TYPE')
        valid_labels = [r[0] for r in cur]
        conn.close()
        return valid_labels

    def seqrecords(self, antibodies, clonal=False, dna=False):
        conn = connect(self.__filename)
        cur = conn.cursor()

        antibodies_ = set(antibodies)

        ab_clause = ' or '.join(['ANTIBODY = ?'] * len(antibodies_))

        equivalencies = set((
            next(cur.execute(
                'select distinct ALT_IDS from ANTIBODY where %s' % ab_clause,
                tuple(antibodies_)
                ))[0]
            or ''
            ).split(',')) - set([''])

        if len(equivalencies):
            antibodies_ |= equivalencies
            ab_clause = ' or '.join(['ANTIBODY = ?'] * len(antibodies_))

        antibodies__ = tuple(sorted(antibodies_))

        stmt = dedent('''\
        select distinct S.NO as NO, S.ID as ID, S.SEQ as SEQ, G.SUBTYPE as SUBTYPE, ? as AB, N.VALUE as VALUE from
        (select SEQUENCE_NO as NO, SEQUENCE_ID as ID, RAW_SEQ as SEQ from SEQUENCE %s group by ID) as S join
        (select SEQUENCE_ID as ID, SUBTYPE from GENO_REPORT group by ID) as G join
        (select SEQUENCE_ID as ID, ANTIBODY as AB, group_concat(TYPE || ':' || VALUE, ',') as VALUE from NEUT where (%s) group by ID) as N
        on N.ID = S.ID and G.ID = S.ID order by S.ID;
        ''' % ('where IS_CLONAL = 1' if clonal else '', ab_clause))
        params = ('+'.join(antibodies__),) + antibodies__
        cur.execute(stmt, params)

        ids = {}
        seqrecords = []
        for row in cur:
            nno, sid, seq, subtype, ab, values = row[:6]
            values_ = {}
            for kv in values.split(','):
                k, v = kv.split(':')
                try:
                    v_ = float(v.strip().lstrip('<>'))
                except ValueError:
                    continue
                if k not in values_:
                    values_[k] = []
                values_[k].append(v_)
            if len(values_) == 0:
                warn("skipping sequence '%s', invalid values '%s'" % (sid, values))
                continue
            dnaseq = Seq(OrfList(seq, include_stops=False)[0], Gapped(generic_nucleotide))
            record = SeqRecord(
                dnaseq if dna else dnaseq.translate(),
                id=sid,
                description=json_dumps({
                    'subtype': subtype,
                    'ab': ab,
                    'values': values_
                    })
            )
            if sid in ids:
                record.id += str(-ids[sid])
                ids[sid] += 1
            else:
                ids[sid] = 1
            seqrecords.append(record)

        conn.close()

        return seqrecords, clonal, antibodies__

    @property
    def subtypes(self):
        conn = connect(self.__filename)
        cur = conn.cursor()
        cur.execute('''select distinct SUBTYPE from GENO_REPORT''')
        valid_subtypes = [r[0] for r in cur if r[0].strip() != '']
        conn.close()
        return valid_subtypes


class MonogramData:
    __no_header_msg = 'input data is not a valid Monogram dataset (no column headers)'

    def __init__(self, fastafile, csvfile):
        self.__fastafile = fastafile
        self.__csvfile = csvfile

    @property
    def basename_root(self):
        return splitext(basename(self.__csvfile))[0]

    @property
    def csvfile(self):
        return self.__csvfile

    @property
    def fastafile(self):
        return self.__fastafile

    @property
    def antibodies(self):
        antibodies = []
        with open(self.__csvfile) as fh:
            sample = fh.read(2048)
            sniffer = csv_sniffer()
            dialect = sniffer.sniff(sample)
            if not sniffer.has_header(sample):
                raise ValueError(MonogramData.__no_header_msg)
            fh.seek(0)
            reader = csv_reader(fh, dialect)
            # grab everything after the accession column in the header row
            for row in reader:
                antibodies.extend(r.strip() for r in row[1:])
                break
        return antibodies

    @property
    def labels(self):
        return []

    def seqrecords(self, antibodies, clonal=False, dna=False):
        if clonal:
            raise ValueError('clonal property is not available with Monogram datasets')
        if dna:
            raise ValueError('dna sequences are not available with Monogram datasets')

        seqrecords = []
        with open(self.__fastafile) as fh:
            seqrecords = [r for r in SeqIO.parse(fh, 'fasta')]

        underdash = re_compile(r'[_-](\d+)$')
        for r in seqrecords:
            r.id = underdash.sub(r'_\1', r.id)

        ic50s = dict((r.id, []) for r in seqrecords)

        with open(self.__csvfile) as fh:
            sniffer = csv_sniffer()
            dialect = sniffer.sniff(fh.read(8192))
            fh.seek(0)
            if not sniffer.has_header(fh.read(8192)):
                raise ValueError(MonogramData.__no_header_msg)
            fh.seek(0)
            reader = csv_reader(fh, dialect)
            columns = None
            for i, row in enumerate(reader):
                if columns is None:
                    columns = dict((v.strip(), j) for j, v in enumerate(row))
                    missing = set(antibodies) - set(columns.keys())
                    if len(missing):
                        raise ValueError("antibodies ('%s') not found!" % "', '".join(missing))
                else:
                    acc = underdash.sub(r'_\1', row[0])
                    try:
                        if acc in ic50s:
                            cln_ic50s = [str(min(float(row[columns[ab]].strip().lstrip('<>')), 25.))
                                         for ab in antibodies
                                         if ab in columns and columns[ab] < len(row)]
                            ic50s[acc].extend(cln_ic50s)
                    except:
                        pass

        drop = []
        for i, r in enumerate(seqrecords):
            if r.id not in ic50s or len(ic50s[r.id]) == 0:
                drop.append(i)
                warn("skipping sequence '%s', VALUE not found" % r.id)
            else:
                r.description = '|'.join(('', '+'.join(antibodies), '', ','.join(ic50s[r.id])))

        for i in sorted(drop, reverse=True):
            del seqrecords[i]

        return seqrecords, clonal, antibodies

    @property
    def subtypes(self):
        return []
