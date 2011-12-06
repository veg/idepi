
import logging, re

from collections import namedtuple
from io import StringIO
from math import copysign, sqrt
from operator import itemgetter
from os import close, remove, rename
from os.path import exists
from shutil import copyfile
from sys import stderr, stdout
from tempfile import mkstemp
from unicodedata import combining

from Bio import AlignIO, SeqIO
from Bio.Alphabet import Gapped, generic_nucleotide, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mrmr import MRMR_LOGGER

from ._hmmer import Hmmer
from ._logging import IDEPI_LOGGER
from ._normalvalue import NormalValue
from ._simulation import Simulation


__all__ = [
    'cv_results_to_output',
    'make_output_meta',
    'crude_sto_read',
    'generate_alignment_from_seqrecords',
    'generate_alignment',
    'pretty_fmt_results',
    'pretty_fmt_stats',
    'pretty_fmt_weights',
    'extract_feature_weights'
]


# refseq_offs has for keys 0-indexed positions into the trimmed refseq
# and for values the number of trimmed positions occuring immediately
# before the key-index. This so the colfilter can properly name the
# positions according to the full reference sequence 
def refseq_off(refseq):
    off = 0
    offs = {}
    i = 0
    for pos in refseq:
        if 'a' < pos and pos < 'z':
            off += 1
        elif 'A' < pos and pos < 'Z':
            if off > 0:
                offs[i] = off
                off = 0
            i += 1
        else:
            continue
    return offs


def crude_sto_read(filename, ref_id_func=None, dna=False):
    Fake = namedtuple('FakeSeqRecord', ['id'])
    alph = Gapped(generic_nucleotide if dna else generic_protein)
    refseq = None
    with open(filename) as fh, StringIO() as tmp:
        notrel = re.compile(r'^(?:#|#=|//)')
        trim = re.compile(r'[^-A-Z]')
        for line in (line.strip() for line in fh):
            if notrel.match(line) or line == '':
                # don't print the GR lines, they won't match anymore after trimming
                if not line.startswith('#=GR'):
                    print(line, file=tmp)
            else:
                padlen = 1
                for i in range(line.find(' ') + 1, len(line)):
                    if line != ' ':
                        break
                    padlen += 1
                pad = ' ' * padlen
                acc, seq = line.split()
                if ref_id_func is not None and ref_id_func(Fake(acc)):
                    refseq = seq
                seq = trim.sub('', seq)
                print(pad.join((acc, seq)), file=tmp)
        tmp.seek(0)
        alignment = AlignIO.read(tmp, 'stockholm', alphabet=alph)

    if refseq is None and ref_id_func is not None:
        raise RuntimeError('Unable to find the reference sequence to compute an offset!')

    offs = None if ref_id_func is None else refseq_off(refseq)
    return alignment, offs


def generate_alignment_from_seqrecords(seq_records, my_basename, opts):
    fd, ab_fasta_filename = mkstemp(); close(fd)
    fd, hmm_filename = mkstemp(); close(fd)
    fd, sto_filename = mkstemp(); close(fd)
    finished = False

    log = logging.getLogger(IDEPI_LOGGER)

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
        except TypeError as e:
            print(seq_records, file=stderr)
            raise e

        # close the handles in this order because BioPython wiki does so
        fafh.close()
        hxb2fh.close()

        # make the tempfiles for the alignments, and close them out after we use them
        SeqIO.convert(opts.REFSEQ_FASTA, 'fasta', sto_filename, 'stockholm')

        hmmer = Hmmer(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)

        numseqs = len(seq_records)

        log.debug('beginning alignment')

        for i in range(opts.HMMER_ITER):
            log.debug('aligning %d sequences (%d of %d)' % (numseqs, i+1, opts.HMMER_ITER))
            hmmer.build(hmm_filename, sto_filename)
            hmmer.align(hmm_filename, ab_fasta_filename, output=sto_filename, alphabet=Hmmer.DNA if opts.DNA else Hmmer.AMINO, outformat=Hmmer.PFAM)

        # rename the final alignment to its destination
        finished = True

    finally:
        # cleanup these files
        if finished:
            copyfile(sto_filename, my_basename + '.sto')
            copyfile(hmm_filename, my_basename + '.hmm')
            log.debug('finished alignment, output moved to %s.sto' % my_basename)
        remove(sto_filename)
        remove(ab_fasta_filename)
        remove(hmm_filename)

def generate_alignment(seqrecords, my_basename, ref_id_func, opts):
    sto_filename = my_basename + '.sto'
    hmm_filename = my_basename + '.hmm'

    if hasattr(opts, 'SIM') and opts.SIM == Simulation.DUMB:
        # we're assuming pre-aligned because they're all generated from the same refseq
        with open(hmm_filename, 'w') as fh:
            SeqIO.write(seqrecords, fh, 'stockholm')
    elif not exists(sto_filename):
        generate_alignment_from_seqrecords(
            seqrecords,
            my_basename,
            opts
        )

    if not exists(hmm_filename):
        hmmer = Hmmer(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)
        hmmer.build(hmm_filename, sto_filename)

    return crude_sto_read(sto_filename, ref_id_func, opts.DNA)


def cv_results_to_output(results, colnames, meta=None):

    statsdict = results.stats.todict()

    # remove minstat 'cause we don't want it here..
    if 'Minstat' in statsdict:
        del statsdict['Minstat']

    featureweights = {}

    for i in range(len(results.extra)):
        assert(len(results.extra[i]['features']) >= len(results.extra[i]['weights']))
        for j in range(len(results.extra[i]['features'])):
            v = results.extra[i]['weights'][j] if j < len(results.extra[i]['weights']) else 0.
            k = results.extra[i]['features'][j]
            if k not in featureweights:
                featureweights[k] = []
            featureweights[k].append(int(copysign(1, v)))

    weightsdict = {}
    for idx, weights in featureweights.items():
        val = NormalValue(int, weights)
        weightsdict[colnames[idx]] = val

    log = logging.getLogger(MRMR_LOGGER)
    log.debug('mrmr index to name map: {%s}' % ', '.join(
        "%d: '%s'" % (
            idx, colnames[idx]
        ) for idx in sorted(featureweights.keys())
    ))

    ret = {}

    if meta is not None:
        ret['meta'] = meta

    ret['statistics'] = dict((k, { 'mean': v.mu, 'std': sqrt(v.sigma) }) for k, v in statsdict.items())
    ret['weights'] = [{ 'position': k, 'value': { 'mean': v.mu, 'std': sqrt(v.sigma), 'N': len(v) } } for k, v in sorted(
        weightsdict.items(),
        key=lambda x: int(re.sub(r'[a-zA-Z\[\]]+', '', x[0]))
    )]

    return ret


def make_output_meta(opts, N, target, antibody, forward_select=None):
    cutoff = opts.IC50LT if target == 'lt' else opts.IC50GT
    return {
        'sequences': N,
        'features': opts.NUM_FEATURES if forward_select is None else forward_select,
        'discriminator': { 'orientation': target, 'cutoff': cutoff },
        'antibody': antibody,
        'folds': opts.CV_FOLDS
    }


def pretty_fmt_stats(stats, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix

    buf += '"statistics": {\n'

    stat_prefixes = {}
    for k in stats.keys():
        stat_prefixes[k] = sum([1 for c in k if combining(c) == 0])

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


def pretty_fmt_weights(weights, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix + '"weights": [\n'

    if len(weights) > 0:
        name_len = max(len(v['position']) for v in weights) + 3
        mean_len = max(len('% .6f' % v['value']['mean']) for v in weights)
        std_len = max(len('%.6f' % v['value']['std']) for v in weights)
        N_len = max(len('%d' % v['value']['N']) for v in weights)
        fmt = '{ "mean": %%%d.6f, "std": %%%d.6f, "N": %%%dd }' % (mean_len, std_len, N_len)
        output = (prefix + '  { "position": %-*s "value": %s }' % (
            name_len, '"%s",' % v['position'],
            fmt % (
                v['value']['mean'],
                v['value']['std'],
                v['value']['N']
            ),
        ) for v in sorted(weights, key=lambda x: int(re.sub(r'[a-zA-Z\[\]]+', '', x['position']))))

    return buf + ',\n'.join(output) + '\n' + prefix + ']'


def pretty_fmt_meta(meta, ident=0):
    prefix = ' ' * 2 * ident

    buf = prefix + '"meta": {\n'

    name_len = max(len(k) for k in meta.keys()) + 3
    output = (prefix + '  %-*s %s' % (
        name_len,
        '"%s":' % k,
        '"%s"' % v if isinstance(v, str) else \
        ' { %s }' % ', '.join(
            ('"%s": %s' % (
                k,
                '%s' % v if isinstance(v, (float, int)) else '"%s"' % v
            ) for k, v in v.items())
        ) if isinstance(v, dict) else \
        ' %s' % str(v)
    ) for k, v in sorted(meta.items(), key=itemgetter(0)))

    return buf + ',\n'.join(output) + '\n' + prefix + '}'


def pretty_fmt_results(results):
    ret = '{\n'
    ret += pretty_fmt_meta(results['meta'], 1) + ',\n' if 'meta' in results else ''
    ret += pretty_fmt_stats(results['statistics'], 1) + ',\n'
    ret += pretty_fmt_weights(results['weights'], 1) + '\n}'
    return ret

def extract_feature_weights(instance):
    return {
        'features': instance.features(),
        'weights': instance.classifier.weights()
    }
