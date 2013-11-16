from __future__ import print_function
import HTSeq
from collections import defaultdict
import sys
from dregs import binf

def extract_exons(fasta_fname, gff_fname):
  sequences = HTSeq.FastaReader(fasta_fname)
  # end_included=True as (exon.end - exon.start) % 3 = 2.
  gff = HTSeq.GFF_Reader(gff_fname, end_included=True)

  features = defaultdict(lambda: defaultdict(list))
  for feat in gff:
    features[feat.name][feat.type].append(feat)

  for kog, feats in features.items():
    exons = feats['Exon']
    exons = sorted(exons, key=lambda e: e.iv.start)
    seq = ''.join([str(sequences[exon.iv]) for exon in exons])
    binf.write_fasta_seq(sys.stdout, kog, seq)

def main():
  if len(sys.argv[1:]) != 2:
    print('Usage: %s [FASTA file] [GFF file]' % sys.argv[0], file=sys.stderr)
    sys.exit(1)
  extract_exons(sys.argv[1], sys.argv[2])

main()
