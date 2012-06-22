
from __future__ import division, print_function

from collections import namedtuple
from json import dumps as json_dumps, loads as json_loads
from logging import getLogger
from math import copysign, sqrt
from operator import itemgetter
from os import close, remove, rename
from os.path import exists
from re import compile as re_compile, sub as re_sub
from shutil import copyfile
from six import u
from sys import stderr, stdout
from tempfile import mkstemp
from unicodedata import combining
from warnings import warn

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import Gapped, generic_nucleotide, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from mrmr import MRMR_LOGGER

from ._hmmer import Hmmer
from ._logging import IDEPI_LOGGER
from ._normalvalue import NormalValue
from ._simulation import Simulation
from ._util import BASE_ALPH, base_10_to_n, base_26_to_alph


__all__ = [
    'cv_results_to_output',
    'make_output_meta',
    'crude_sto_read',
    'generate_alignment_from_seqrecords',
    'generate_alignment',
    'pretty_fmt_meta',
    'pretty_fmt_results',
    'pretty_fmt_stats',
    'pretty_fmt_weights',
    'extract_feature_weights_similar',
    'extract_feature_weights',
    'column_labels',
    'fasta_json_desc',
    'fix_hxb2_seq'
]


try:
    unicode is not None
except NameError:
    unicode = str


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


def crude_sto_read(filename, ref_id_func=None, dna=False, description=False):
    Fake = namedtuple('FakeSeqRecord', ['id'])
    alph = Gapped(generic_nucleotide if dna else generic_protein)
    refseq = None
    msa = MultipleSeqAlignment([], alphabet=alph)
    with open(filename) as fh:
        notrel = re_compile(r'^(?://|#(?!=GS))')
        isdesc = re_compile(r'^#=GS')
        trim = re_compile(r'[^-A-Z]+')
        descs = {}
        for line in fh:
            line = line.strip()
            if line == '' or notrel.match(line):
                continue
            elif isdesc.match(line):
                if description:
                    elems = line.split(None, 3)
                    if len(elems) > 3 and elems[2] == 'DE':
                        acc, desc = elems[1], elems[3]
                        if acc in descs:
                            warn("duplicate sequence name '%s' detected! The stockholm specification doesn't allow this!" % acc)
                        descs[acc] = desc
                else:
                     continue
            else:
                try:
                    acc, seq = line.split(None, 1)
                except ValueError:
                    warn("skipping line '%s', doesn't seem to contain a sequence" % line)
                    continue
                if ref_id_func is not None and ref_id_func(Fake(acc)):
                    refseq = seq
                try:
                    seq = trim.sub('', seq)
                    desc = descs[acc] if acc in descs else acc
                    msa.append(SeqRecord(Seq(seq), id=acc, description=desc))
                except ValueError:
                    warn("skipping sequence '%s', it doesn't match the length of the MSA (%d vs %d)" % (acc, len(seq), msa.get_alignment_length()))

    if refseq is None and ref_id_func is not None:
        raise RuntimeError('Unable to find the reference sequence to compute an offset!')

    offs = None if ref_id_func is None else refseq_off(refseq)
    return msa, offs


def generate_alignment_from_seqrecords(seq_records, my_basename, opts):
    fd, ab_fasta_filename = mkstemp(); close(fd)
    fd, hmm_filename = mkstemp(); close(fd)
    fd, sto_filename = mkstemp(); close(fd)
    finished = False

    log = getLogger(IDEPI_LOGGER)

    try:
        # get the FASTA format file so we can HMMER it
        fafh = open(ab_fasta_filename, 'w')

        # insert the reference sequence
        seq_records.append(opts.REFSEQ)

        # type errors here?
        try:
            SeqIO.write(seq_records, fafh, 'fasta')
        except TypeError as e:
            print(seq_records, file=stderr)
            raise e

        # close the handles in this order because BioPython wiki does so
        fafh.close()

        # make the tempfiles for the alignments, and close them out after we use them
        with open(sto_filename, 'w') as fh:
            SeqIO.write([opts.REFSEQ], fh, 'stockholm')

        hmmer = Hmmer(opts.HMMER_ALIGN_BIN, opts.HMMER_BUILD_BIN)

        numseqs = len(seq_records)

        log.debug('beginning alignment')

        for i in range(opts.HMMER_ITER):
            log.debug('aligning %d sequences (%d of %d)' % (numseqs, i+1, opts.HMMER_ITER))
            hmmer.build(hmm_filename, sto_filename)
            hmmer.align(
                hmm_filename,
                ab_fasta_filename,
                output=sto_filename,
                alphabet=Hmmer.DNA if opts.DNA else Hmmer.AMINO,
                outformat=Hmmer.PFAM
            )

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

    # we store the antibody information in the description, so grab it
    return crude_sto_read(sto_filename, ref_id_func, opts.DNA, description=True)


def cv_results_to_output(results, colnames, meta=None, similar=True):

    statsdict = results.stats.todict()

    # remove minstat 'cause we don't want it here..
    if 'Minstat' in statsdict:
        del statsdict['Minstat']

    idxnames = {}
    weightsdict = {}

    for i in range(len(results.extra)):
        # weights - 1 to account for bias term we added in LinearSvm.__bias
        assert(len(results.extra[i]['features']) >= len(results.extra[i]['weights']) - 1)
        # we only ever go up to the # of features selected by mRMR,
        # so the last weight (of the bias term) is ignored here, intentionally
        for j in range(len(results.extra[i]['features'])):
            v = results.extra[i]['weights'][j] if j < len(results.extra[i]['weights']) else 0.
            k = results.extra[i]['features'][j]
            r = set([idx for idx, _ in results.extra[i]['similar'][k]]) if similar else None
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

    ret = {}

    if meta is not None:
        ret['meta'] = meta

    numeric = re_compile(r'[^0-9]+')

    ret['statistics'] = dict((k, { 'mean': v.mu, 'std': sqrt(v.sigma) }) for k, v in statsdict.items())
    ret['weights'] = [{ 'position': k, 'value': { 'mean': v[0].mu, 'std': sqrt(v[0].sigma), 'N': len(v[0]) } } for k, v in sorted(
        weightsdict.items(),
        key=lambda x: int(numeric.sub('', x[0]))
    )]

    if similar:
        for i, kv in enumerate(sorted(weightsdict.items(), key=lambda x: int(numeric.sub('', x[0])))):
            _, v = kv
            ret['weights'][i]['similar'] = [colnames[j] for j in v[1]]

    return ret


def make_output_meta(opts, N, balance, target, antibody, forward_select=None):
    cutoff = opts.IC50LT if target in ('le', 'lt') else opts.IC50GT
    return {
        'sequences': N,
        'balance': balance,
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


def pretty_fmt_weights(weights, ident=0, similar=True):
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
                '"%s"' % v if isinstance(v, str) else
                '%.6g' % v if isinstance(v, float) else
                '%s' % str(v)
            ) for k, v in v.items())
        ) if isinstance(v, dict) else \
        ' %.6g' % v if isinstance(v, float) else \
        ' %s' % str(v)
    ) for k, v in sorted(meta.items(), key=itemgetter(0)))

    return buf + ',\n'.join(output) + '\n' + prefix + '}'


def pretty_fmt_results(results, similar=True):
    ret = '{\n'
    ret += pretty_fmt_meta(results['meta'], 1) + ',\n' if 'meta' in results else ''
    ret += pretty_fmt_stats(results['statistics'], 1) + ',\n'
    ret += pretty_fmt_weights(results['weights'], 1, similar) + '\n}'
    return ret


def extract_feature_weights_similar(instance, similar=True):
    ret = {
        'features': instance.features(),
        'weights':  instance.classifier.weights()
    }
    if similar:
        ret['similar'] = instance.selector.related()
    return ret


def extract_feature_weights(instance):
    return extract_feature_weights_similar(instance, False)


def column_labels(refseq, refseq_offs={}):
    if isinstance(refseq, str):
        ref = refseq
    elif isinstance(refseq, SeqRecord):
        ref = str(refseq.seq)
    colnames = []
    colnum = 0
    insert = 0
    for i, p in enumerate(ref):
        if i in refseq_offs:
            colnum += refseq_offs[i]
        if ref[i] not in '._-':
            colnum += 1
            insert = 0
        else:
            insert += 1
        colname = '%s%d%s' % (p if insert == 0 else '', colnum, base_26_to_alph(base_10_to_n(insert, BASE_ALPH)))
        colnames.append(colname)
    return colnames


def fasta_json_desc(seqrecord):
    try:
        return json_loads(seqrecord.description.strip(seqrecord.id).strip())
    except ValueError:
        return {}


def fix_hxb2_seq(opts):
    from BioExt import hxb2
    if opts.DNA == False and str(opts.REFSEQ.seq) == str(hxb2.env.load().seq):
        try:
            opts.REFSEQ.seq = opts.REFSEQ.seq.translate()
            data = fasta_json_desc(opts.REFSEQ)
            if isinstance(data, dict) and 'loops' in data:
                for k in data['loops'].keys():
                    for i, v in enumerate(data['loops'][k]):
                        data['loops'][k][i] = int(v // 3)
                opts.REFSEQ.description = ' '.join((opts.REFSEQ.id, json_dumps(data, separators=(',', ':'))))
        except ValueError:
            pass # we're already a protein
