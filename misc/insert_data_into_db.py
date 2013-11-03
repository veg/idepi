

import sqlite
import sys

from argparse import ArgumentParser

from idepi.argument import PathType
from idepi.datasource import DataSource


def main(args):

    parser = ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--csv',    type=PathType, dest='source', nargs=2, metavar=('FASTA', 'CSV'))
    group.add_argument('--sqlite', type=PathType, dest='source', nargs=1, metavar='SRC_SQLITE3')
    parser.add_argument('--db',    type=PathType, dest='db', nargs=1, metavar='DEST_SQLITE3')

    ns, args = parser.parse_args(args)

    source = DataSource(*ns.source)
    conn = sqlite3.connect(ns.db)
    cur = sqlite3.cursor()

    db_abs = set(cur.execute('select distinct ANTIBODY from ANTIBODY'))
    db_seqs = dict(cur.execute('select distinct RAW_SEQ,SEQUENCE_ID from SEQUENCE'))

    raise RuntimeError('this script is unfinished')

    for ab in source.antibodies:
        records = source.seqrecords(ab)
        pass

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
