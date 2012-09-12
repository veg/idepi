
from __future__ import division, print_function

from csv import reader as csv_reader, Sniffer as csv_sniffer
from logging import getLogger
from os.path import basename, splitext
from re import compile as re_compile
from sqlite3 import OperationalError, connect
from warnings import warn

from Bio import SeqIO
from Bio.Alphabet import Gapped, generic_nucleotide
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from BioExt import OrfList

from ..logging import IDEPI_LOGGER


__all__ = ['DataSource']


_equivalencies = {}
for k, v in [('PG9', 'NAC17'), ('PG16', 'NAC18')]:
    _equivalencies[k] = [v]
    _equivalencies[v] = [k]


def DataSource(s):
    args = s.split(',')
    if len(args) == 1:
        return Sqlite3Db(*args)
    if len(args) == 2:
        return MonogramData(*args)
    else:
        raise ValueError("invalid data specification: '%s'" % s)


class Sqlite3Db(object):

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
        curr = conn.cursor()
        curr.execute('''select distinct ANTIBODY from NEUT''')
        valid_antibodies = [r[0] for r in curr]
        conn.close()
        return valid_antibodies

    def seqrecords(self, antibody, clonal=False, dna=False):
        conn = connect(self.__filename)
        cur = conn.cursor()

        if antibody in _equivalencies:
            antibodies = tuple([antibody, antibody] + _equivalencies[antibody])
            ab_clause = ' OR '.join(['ANTIBODY = ?'] * len(antibodies[1:]))
        else:
            antibodies = (antibody, antibody,)
            ab_clause = 'ANTIBODY = ?'

        try:
            cur.execute('''
            select distinct S.NO as NO, S.ID as ID, S.SEQ as SEQ, G.SUBTYPE as SUBTYPE, ? as AB, N.IC50 as IC50 from
            (select SEQUENCE_NO as NO, SEQUENCE_ID as ID, RAW_SEQ as SEQ from SEQUENCE %s group by ID) as S join
            (select SEQUENCE_ID as ID, SUBTYPE from GENO_REPORT group by ID) as G join
            (select SEQUENCE_ID as ID, ANTIBODY as AB, group_concat(IC50, ',') as IC50 from NEUT where %s group by ID) as N
            on N.ID = S.ID and G.ID = S.ID order by S.ID;
            ''' % ('where IS_CLONAL = 1' if clonal else '', ab_clause), antibodies)
        except OperationalError:
            getLogger(IDEPI_LOGGER).debug('falling back to older database format, CLONAL feature is unsupported')
            clonal = None
            cur.execute('''
            select distinct S.NO as NO, S.ID as ID, S.SEQ as SEQ, G.SUBTYPE as SUBTYPE, ? as AB, N.IC50 as IC50 from
            (select SEQUENCE_NO as NO, ACCESSION_ID as ID, RAW_SEQ as SEQ from SEQUENCE group by ID) as S join
            (select ACCESSION_ID as ID, SUBTYPE from GENO_REPORT group by ID) as G join
            (select ACCESSION_ID as ID, ANTIBODY as AB, group_concat(IC50_STRING, ',') as IC50 from NEUT where %s group by ID) as N
            on N.ID = S.ID and G.ID = S.ID order by S.ID;
            ''' % ab_clause, antibodies)

        ids = {}
        seqrecords = []
        for row in cur:
            nno, sid, seq, subtype, ab, ic50s = row[:6]
            cln_ic50s = []
            for ic50 in ic50s.split(','):
                try:
                    v = min(float(ic50.strip().lstrip('<>')), 25.)
                except ValueError:
                    continue
                cln_ic50s.append(str(v))
            if len(cln_ic50s) == 0:
                warn("skipping sequence '%s', invalid IC50s '%s'" % (sid, ic50s))
                continue
            dnaseq = Seq(OrfList(seq, include_stops=False)[0], Gapped(generic_nucleotide))
            record = SeqRecord(
                dnaseq if dna else dnaseq.translate(),
                id=sid,
                description='|'.join((subtype, ab, ','.join(cln_ic50s)))
            )
            if sid in ids:
                record.id += str(-ids[sid])
                ids[sid] += 1
            else:
                ids[sid] = 1
            seqrecords.append(record)

        conn.close()

        return seqrecords, clonal

    @property
    def subtypes(self):
        conn = connect(self.__filename)
        curr = conn.cursor()
        curr.execute('''select distinct SUBTYPE from GENO_REPORT''')
        valid_subtypes = [r[0] for r in curr if r[0].strip() != '']
        conn.close()
        return valid_subtypes


class MonogramData(object):
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
            sniffer = csv_sniffer()
            dialect = sniffer.sniff(fh.read(8192))
            fh.seek(0)
            if not sniffer.has_header(fh.read(8192)):
                raise ValueError(MonogramData.__no_header_msg)
            fh.seek(0)
            reader = csv_reader(fh, dialect)
            # grab everything after the accession column in the header row
            for row in reader:
                antibodies.extend(r.strip() for r in row[1:])
                break
        return antibodies

    def seqrecords(self, antibody, clonal=False, dna=False):
        if clonal:
            raise ValueError('clonal property is not available with Monogram datasets')
        if dna:
            raise ValueError('dna sequences are not available with Monogram datasets')

        antibodies = [antibody]
        if antibody in _equivalencies:
            antibodies.append(_equivalencies[antibody])

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
                    if antibodies[0] not in columns:
                        raise ValueError("antibody ('%s') not found!" % antibodies[0])
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
                warn("skipping sequence '%s', IC50 not found" % r.id)
            else:
                r.description = '|'.join(('', antibody, ','.join(ic50s[r.id])))

        for i in sorted(drop, reverse=True):
            del seqrecords[i]

        return seqrecords, clonal

    @property
    def subtypes(self):
        return []
