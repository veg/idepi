
import re

from math import copysign, sqrt
from operator import itemgetter
from os import close, remove, rename
from os.path import exists
from shutil import copyfile
from sys import stderr, stdout
from tempfile import mkstemp
from unicodedata import combining

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import Gapped, generic_nucleotide, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from _hmmer import Hmmer
from _normalvalue import NormalValue
from _simulation import Simulation


__all__ = [
    'cv_results_to_output',
    'crude_sto_read',
    'generate_alignment_from_SeqRecords',
    'generate_alignment',
    'pretty_fmt_results',
    'pretty_fmt_stats',
    'pretty_fmt_weights'
]


def crude_sto_read(filename, dna=False):
    with open(filename) as fh:
        notrel = re.compile(r'^(?:#|#=|//)')
        alndict = dict([line.strip().split() for line in fh if not notrel.match(line.strip()) and line.strip() != ''])
    trim = re.compile(r'[^-A-Z]')
    alph = Gapped(generic_nucleotide if dna else generic_protein)
    alignment = MultipleSeqAlignment([
        SeqRecord(
            Seq(trim.sub('', seq), alph),
            acc
        ) for acc, seq in alndict.items()
    ])
    return alignment


def generate_alignment_from_SeqRecords(seq_records, my_basename, opts):
    fd, ab_fasta_filename = mkstemp(); close(fd)
    fd, hmm_filename = mkstemp(); close(fd)
    fd, sto_filename = mkstemp(); close(fd)
    finished = False

    try:
        # get the FASTA format file so we can HMMER it
        fafh = open(ab_fasta_filename, 'w')
        hxb2fh = open(opts.REFSEQ_FASTA, 'rU')

        # grab the HXB2 Reference Sequence
        hxb2_record = SeqIO.parse(hxb2fh, 'fasta')
        seq_records.extend(hxb2_record)

        # type errors here?
        try:
            SeqIO.write(seq_records, fafh, 'fasta')
        except TypeError, e:
            print >> stderr, seq_records
            raise e

        # close the handles in this order because BioPython wiki does so
        fafh.close()
        hxb2fh.close()

        # make the tempfiles for the alignments, and close them out after we use them
        SeqIO.convert(opts.REFSEQ_FASTA, 'fasta', sto_filename, 'stockholm')

        print >> stderr, 'Aligning %d sequences with HMMER:' % (len(seq_records)+1),

        hmmer = Hmmer(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)

        for i in xrange(0, opts.HMMER_ITER):
            print >> stderr, '%d,' % i,
            hmmer.build(hmm_filename, sto_filename)
            hmmer.align(hmm_filename, ab_fasta_filename, output=sto_filename, alphabet=Hmmer.DNA if opts.DNA else Hmmer.AMINO, outformat=Hmmer.PFAM)

        # rename the final alignment to its destination
        finished = True

    finally:
        # cleanup these files
        if finished:
            copyfile(sto_filename, my_basename + '.sto')
            copyfile(hmm_filename, my_basename + '.hmm')
            print >> stderr, 'done, output moved to: %s' % my_basename + '.sto'
        remove(sto_filename)
        remove(ab_fasta_filename)
        remove(hmm_filename)

def generate_alignment(seqrecords, my_basename, opts):
    sto_filename = my_basename + '.sto'
    hmm_filename = my_basename + '.hmm'

    if hasattr(opts, 'SIM') and opts.SIM == Simulation.DUMB:
        # we're assuming pre-aligned because they're all generated from the same refseq
        with open(hmm_filename, 'w') as fh:
            SeqIO.write(seqrecords, fh, 'stockholm')
    elif not exists(sto_filename):
        generate_alignment_from_SeqRecords(
            seqrecords,
            my_basename,
            opts
        )

    if not exists(hmm_filename):
        hmmer = Hmmer(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)
        hmmer.build(hmm_filename, sto_filename)

    alignment = crude_sto_read(sto_filename, opts.DNA)

    return alignment

def cv_results_to_output(results, colnames):

    statsdict = results['stats'].todict()

    # remove minstat 'cause we don't want it here..
    if 'Minstat' in statsdict:
        del statsdict['Minstat']

    featureweights = {}
    for i in xrange(len(results['extra'])):
        assert(len(results['extra'][i]['features']) >= len(results['extra'][i]['weights']))
        for j in xrange(len(results['extra'][i]['features'])):
            v = results['extra'][i]['weights'][j] if j < len(results['extra'][i]['weights']) else 0.
            k = results['extra'][i]['features'][j]
            if k not in featureweights:
                featureweights[k] = []
            featureweights[k].append(int(copysign(1, v)))

    weightsdict = {}
    for idx, weights in featureweights.items():
        val = NormalValue(int, weights)
        weightsdict[colnames[idx]] = val

    ret = {}

    ret['statistics'] = dict([(k, { 'mean': v.mu, 'std': sqrt(v.sigma) }) for k, v in statsdict.items()])
    ret['weights'] = [{ 'position': k, 'value': { 'mean': v.mu, 'std': sqrt(v.sigma), 'N': len(v) } } for k, v in sorted(
        weightsdict.items(),
        key=lambda x: int(re.sub(r'[a-zA-Z\[\]]+', '', x[0]))
    )]

    return ret

def pretty_fmt_stats(stats, ident=0):
    prefix = u' ' * 2 * ident

    buf = prefix

    buf += '"statistics": {\n'

    stat_prefixes = {}
    for k in stats.keys():
        stat_prefixes[k] = sum([1 for c in k if combining(c) == 0])

    stat_len = max(stat_prefixes.values())
    mean_len = max([len('%.6f' % v['mean']) for v in stats.values()])
    std_len = max([len('%.6f' % v['std']) for v in stats.values()])
    fmt = u'{ "mean": %%%d.6f, "std": %%%d.6f }' % (mean_len, std_len)
    output = [prefix + u'  %s%s %s' % (
        u'"%s":' % k,
        u' ' * (stat_len - stat_prefixes[k]),
        fmt % (v['mean'], v['std'])
    ) for k, v in sorted(stats.items(), key=itemgetter(0))]

    return buf + ',\n'.join(output) + '\n' + prefix + '}'

def pretty_fmt_weights(weights, ident=0):
    prefix = u' ' * 2 * ident

    buf = prefix + '"weights": [\n'

    if len(weights) > 0:
        name_len = max([len(v['position']) for v in weights]) + 3
        mean_len = max([len('% .6f' % v['value']['mean']) for v in weights])
        std_len = max([len('%.6f' % v['value']['std']) for v in weights])
        N_len = max([len('%d' % v['value']['N']) for v in weights])
        fmt = u'{ "mean": %%%d.6f, "std": %%%d.6f, "N": %%%dd }' % (mean_len, std_len, N_len)
        output = [prefix + u'  { "position": %-*s "value": %s }' % (
            name_len, u'"%s",' % v['position'],
            fmt % (
                v['value']['mean'],
                v['value']['std'],
                v['value']['N']
            ),
        ) for v in sorted(weights, key=lambda x: int(re.sub(r'[a-zA-Z\[\]]+', '', x['position'])))]

    return buf + ',\n'.join(output) + '\n' + prefix + ']'

def pretty_fmt_results(ret):
    return '{\n' + pretty_fmt_stats(ret['statistics'], 1) + \
           ',\n' + pretty_fmt_weights(ret['weights'], 1) + '\n}'
