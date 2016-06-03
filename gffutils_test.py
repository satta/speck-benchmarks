#!/usr/bin/env python

import sys
import os
import operator
import gffutils
import collections
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

STOPCODONPAT = re.compile('.*[*#+].*')
fails = collections.defaultdict(lambda: collections.defaultdict(int))

def get_protein(fdb, mrna, seqs):
    seq_exons = []
    cdslist = list(fdb.children(mrna, featuretype='CDS'))
    cdslist.sort(key=lambda f: f.start)
    for c in cdslist:
        seq_exons.append(str(seqs[c.seqid][c.start-1:c.end].seq))
    gene_seq = Seq(''.join(seq_exons))
    if mrna.strand == '-':
        gene_seq = gene_seq.reverse_complement()
    return str(gene_seq.translate())

def has_transcript(fdb, gene, seqs):
    if sum(1 for i in
           fdb.children(gene, featuretype=['mRNA','tRNA','snRNA',
                                           'rRNA','snoRNA','ncRNA'])) == 0:
        fails['gene']['has_transcript'] += 1

def children_in_coords(fdb, gene, seqs):
        for c in fdb.children(gene):
            if c.start < gene.start or c.end > gene.end:
                fails['gene']['children_in_coords'] += 1

def children_consistent_strands(fdb, gene, seqs):
    for c in fdb.children(gene):
        if c.strand != gene.strand:
            fails['gene']['children_consistent_strands'] += 1

def not_suspiciously_short(fdb, gene, seqs):
    if gene.end - gene.start + 1 <= 30:
        fails['gene']['not_suspiciously_short'] += 1

def n_content(fdb, mrna, seqs):
    mrna_seq = seqs[mrna.seqid][mrna.start-1:mrna.end].seq.lower()
    if mrna_seq.count('n')/float(mrna.end - mrna.start + 1) >= 0.5:
        fails['mRNA']['n_content'] += 1

def has_CDS(fdb, mrna, seqs):
    if sum(1 for i in fdb.children(mrna, featuretype='CDS')) == 0:
        fails['mRNA']['has_CDS'] += 1

def has_only_CDS_and_UTR_children(fdb, mrna, seqs):
    num_children = sum(1 for i in fdb.children(mrna))
    num_cds_utr = sum(1 for i in fdb.children(mrna,
                        featuretype=['CDS','five_prime_UTR','three_prime_UTR']))
    if num_children != num_cds_utr:
        fails['mRNA']['has_only_CDS_and_UTR_children'] += 1

def minimum_length(fdb, mrna, seqs):
    if mrna.end - mrna.start + 1 < 3:
        fails['mrna']['minimum_length'] += 1

def internal_stop_codon(fdb, mrna, seqs):
    prot_seq = get_protein(fdb, mrna, seqs)
    if STOPCODONPAT.match(prot_seq[:-1]):
        fails['mrna']['internal_stop_codon'] += 1

def has_no_children(fdb, cds, seqs):
    if sum(1 for i in fdb.children(cds)) > 0:
        fails['cds']['has_no_children'] += 1

checks = {}
checks = { 'gene' : [has_transcript,
                     children_in_coords,
                     children_consistent_strands,
                     not_suspiciously_short],
           'mRNA' : [n_content,
                     has_CDS,
                     has_only_CDS_and_UTR_children,
                     minimum_length,
                     internal_stop_codon],
           'CDS'  : [has_no_children] }

def main(gff_file, ref_file):
    seqs = SeqIO.to_dict(SeqIO.parse(ref_file, "fasta"))

    fdb = gffutils.FeatureDB(gff_file, keep_order=True)
    #print get_protein(fdb, fdb['LmjF.12.0867:mRNA'], seqs)
    #sys.exit()

    for f in fdb.all_features():
        if checks.has_key(f.featuretype):
            for check in checks[f.featuretype]:
                check(fdb, f, seqs)
    for x in fails:
        print(x)
        for y in fails[x]:
            print (y, fails[x][y])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit()
    main(*sys.argv[1:])
