
from __future__ import division, print_function

import json
from gzip import GzipFile
from io import StringIO
from Bio import SeqIO


__all__ = ['PhyloGzFile']


class PhyloGzFile:

    def __init__(self):
        pass

    def read(filename):
        tree = None
        alignment = None
        colnames = None
        xdata = None

        with StringIO() as tmp:
            with GzipFile(filename) as fh:
                tmp.write(fh.read().decode('utf-8'))
            tmp.seek(0)
            input_dict = json.load(tmp)

        for part in ('tree', 'alignment', 'colnames', 'xdata'):
            if part not in input_dict:
                raise RuntimeError('%s datum not found in PhyloGzFile %s' % (part, filename))

        tree = input_dict['tree']

        with StringIO(input_dict['alignment']) as tmp:
            alignment = [r for r in SeqIO.parse(tmp, 'fasta')]

        colnames = input_dict['colnames']

        xdata = input_dict['xdata']

        assert (tree is not None), 'Tree not found!'
        assert (alignment is not None), 'Alignment not found!'
        assert (colnames is not None), 'Column labels not found!'
        assert (xdata is not None), 'Extra data not found!'

        return tree, alignment, colnames, xdata


    def write(filename, tree, alignment, colnames, xdata={}):
        if len(filename) < 4 or filename[-4:].lower() != '.pgz':
            filename += '.pgz'

        with GzipFile(filename, 'wb') as fh:
            output_dict = { 'tree': tree, 'colnames': colnames, 'xdata': xdata }

            with StringIO() as tmp:
                SeqIO.write(alignment, tmp, 'fasta')
                output_dict['alignment'] = tmp.getvalue().strip()

            json.dump(output_dict, fh, separators=(',', ':'))
