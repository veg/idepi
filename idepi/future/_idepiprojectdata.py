
from __future__ import division, print_function

import bz2
import json
import sqlite3

from datetime import datetime
from os import close
from os.path import exists
from six import StringIO
from tempfile import mkstemp

import numpy as np


__all__ = ['IdepiProjectData']


class IdepiProjectData:
    __HMM = 'hmm'
    __ALIGNMENT = 'alignment'
    __COLUMNS = 'columns'
    __DISCRETE_DATA = 'discretedata'
    __PHYLO_DATA = 'phylodata'
    __META = 'meta'
    __RUNS = 'runs'

    def __create(self):
        if not exists(self.__filename):
            conn = sqlite3.connect(self.__filename)
            cur = conn.cursor()
            cur.execute('create table %s (id INTEGER PRIMARY KEY ASC, latest_discrete TEXT UNIQUE, latest_phylo TEXT UNIQUE, antibody TEXT UNIQUE' % IdepiProjectData.__META)
            cur.execute('create table %s (id INTEGER PRIMARY KEY ASC, key TEXT UNIQUE, date TEXT, result BLOB)' % IdepiProjectData.__RUNS)
            for table in (
                    IdepiProjectData.__HMM,
                    IdepiProjectData.__ALIGNMENT,
                    IdepiProjectData.__COLUMNS,
                    IdepiProjectData.__DISCRETE_DATA,
            ):
                cur.execute('create table %s (id INTEGER PRIMARY KEY ASC, key TEXT UNIQUE, data BLOB)' % table)
                cur.execute('create unique index idx_%s on %s (key ASC)' % (table, table))

    def __init__(self, filename):
        self.__filename = filename

    def save(self, key, hmmfile, stofile, colnames, x, results):
        datatable = IdepiProjectData.__DISCRETE_DATA
        if x.dtype == float:
            datatable = IdepiProjectData.__PHYLO_DATA
        conn = sqlite3.connect(self.__filename)
        cur = conn.cursor()
        cur.execute('select key from %s where key = ?' % datatable, (key,))
        if not cur.rowcount:
            self._save_data(key, hmmfile, stofile, colnames, x, cur)
        self._save_run(key, results, cur, x.dtype == float)
        conn.close()

    def _save_data(self, key, hmmfile, stofile, colnames, x, cur=None):
        with open(hmmfile) as fh:
            hmmdata = fh.read()

        with open(stofile) as fh:
            stodata = fh.read()

        phylo = True if x.dtype == float else False

        coldata = json.dumps(colnames)
        xfile = StringIO()
        np.save(xfile, x)
        xdata = bz2.compress(xfile.getvalue(), 9)
        xfile.close()

        closable = False
        if cur is None:
            conn = sqlite3.connect(self.__filename)
            cur = conn.cursor()
            closable = True

        IdepiProjectData.__insert(cur, IdepiProjectData.__HMM, key, hmmdata)
        IdepiProjectData.__insert(cur, IdepiProjectData.__ALIGNMENT, key, stodata)
        IdepiProjectData.__insert(cur, IdepiProjectData.__COLUMNS, key, coldata)
        if phylo:
            IdepiProjectData.__insert(cur, IdepiProjectData.__PHYLO_DATA, key, xdata)
        else:
            IdepiProjectData.__insert(cur, IdepiProjectData.__DISCRETE_DATA, key, xdata)

        if closable:
            conn.close()

    def _save_run(self, key, results, cur=None, phylo=False):
        closable = False
        if cur is None:
            conn = sqlite3.connect(self.__filename)
            cur = conn.cursor()
            closable = True

        cur.execute('insert into %s values (?, ?, ?)' % IdepiProjectData.__RUNS, (key, datetime.utcnow().isoformat(), json.dumps(results)))
        cur.execute('update meta %s = ? where id = 0' % 'latest_phylo' if phylo else 'latest_discrete', (key,))

        if closable:
            conn.close()

    @staticmethod
    def __insert(cur, table, key, value):
        cur.execute('insert into %s values (?, ?)' % table, (key, value))

    @staticmethod
    def __select(cur, table, key):
        cur.execute('select data from %s where key = ?' % table, key)
        if not cur.rowcount:
            raise ValueError('Something broke while trying to read the project file')
        return cur.fetchone()

    def load(self, phylo=False):
        datatable = IdepiProjectData.__DISCRETE_DATA
        if phylo:
           datatable = IdepiProjectData.__PHYLO_DATA
        conn = sqlite3.connect(self.__filename)
        cur = conn.cursor()
        cur.execute('select %s from %s' % ('latest_phylo' if phylo else 'latest_discrete', IdepiProjectData.__META))
        if not cur.rowcount:
            return None
        key = cur.fetchone()
        fd, hmmfile = mkstemp(); close(fd)
        fd, stofile = mkstemp(); close(fd)
        hmmdata = IdepiProjectData.__select(cur, IdepiProjectData.__HMM, key)
        stodata = IdepiProjectData.__select(cur, IdepiProjectData.__ALIGNMENT, key)
        colnames = json.loads(IdepiProjectData.__select(cur, IdepiProjectData.__COLUMNS, key))
        xfile = StringIO()
        xfile.write(bz2.decompress(IdepiProjectData.__select(cur, datatable, key)))
        xfile.seek(0)
        x = np.load(xfile)
        xfile.close()
        conn.close()
        return key, hmmdata, stodata, colnames, x
