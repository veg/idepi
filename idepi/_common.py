
from math import copysign, sqrt
from os import close, remove, rename
from os.path import exists
from re import sub
from sys import stderr
from tempfile import mkstemp

from Bio import AlignIO, SeqIO

from _hmmer import Hmmer
from _normalvalue import NormalValue
from _simulation import Simulation


__all__ = [
    'cv_results_to_output',
    'generate_alignment_from_SeqRecords',
    'generate_alignment'
]


def generate_alignment_from_SeqRecords(seq_records, filename, opts):
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
        print >> stderr, 'done, output moved to: %s' % filename
        finished = True

    finally:
        # cleanup these files
        if finished:
            rename(sto_filename, filename)
        else:
            remove(sto_filename)
        remove(ab_fasta_filename)
        remove(hmm_filename)

def generate_alignment(seqrecords, filename, opts):
    if hasattr(opts, 'SIM') and opts.SIM == Simulation.DUMB:
        # we're assuming pre-aligned because they're all generated from the same refseq
        fh = open(filename, 'w')
        SeqIO.write(seqrecords, fh, 'stockholm')
        fh.close()
    elif not exists(filename):
        generate_alignment_from_SeqRecords(
            seqrecords,
            opts,
            filename
        )

    with open(filename) as fh:
        alignment = AlignIO.read(fh, 'stockholm')

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

    ret['statistics'] = dict([(k.lower(), { 'mean': v.mu, 'std': sqrt(v.sigma) }) for k, v in statsdict.items()])
    ret['weights'] = [{ 'position': k, 'value': { 'mean': v.mu, 'std': sqrt(v.sigma), 'N': len(v) } } for k, v in sorted(
        weightsdict.items(),
        key=lambda x: int(sub(r'[a-zA-Z\[\]]+', '', x[0]))
    )]

    return ret

