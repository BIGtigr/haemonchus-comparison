#!/usr/bin/env python2
'''
Filter isoforms, retaining only the longest amino acid sequence for each gene.

Usage: filter_isoforms.py <input FASTA file> <delimiter>

For each FASTA header, the script will split on whitespace, then take the first
field as the sequence ID. It will set the gene ID to be everything preceding
the *last* occurrence of <delimiter> in the sequence ID. For each distinct gene
ID, it will then output the longest corresponding sequence from amongst the set
sharing the same gene ID.
'''

import sys
from collections import defaultdict
from dregs import binf

def determine_gene_and_isoform_ids(seq_id, delimiter):
  gene_id, isoform_id = seq_id.rsplit(delimiter, 1)
  return (gene_id, isoform_id)

def extract_ids(header, delimiter):
  tokens = header.split()
  seq_id = tokens[0]
  gene_id, isoform_id = determine_gene_and_isoform_ids(seq_id, delimiter)
  return (seq_id, gene_id, isoform_id)

def calculate_sequence_lengths(in_file, delimiter):
  seq_lengths = defaultdict(lambda: defaultdict(int))

  for header, seq in binf.parse_fasta(in_file):
    seq_id, gene_id, isoform_id = extract_ids(header, delimiter)
    seq_lengths[gene_id][isoform_id] = len(seq)

  return seq_lengths

def find_max_isoform_length_for_gene(gene_id, seq_lengths):
  return max(seq_lengths[gene_id].values())

def print_longest_isoforms_for_each_gene(in_file, delimiter, seq_lengths):
  for header, seq in binf.parse_fasta(in_file):
    seq_id, gene_id, isoform_id = extract_ids(header, delimiter)
    if gene_id not in seq_lengths:
      continue

    max_isoform_len = find_max_isoform_length_for_gene(gene_id, seq_lengths)

    if len(seq) == max_isoform_len:
      binf.write_fasta_seq(sys.stdout, header, seq)
      # To mark a given gene as already having had its longest isoform printed,
      # remove the gene from seq_lengths.
      del seq_lengths[gene_id]

def main():
  # Can't use stdin, as must read *twice* -- first to calculate length of each sequence, and then to output longest sequence for each.
  input_fasta_filename = sys.argv[1]

  # Delimiter: character after whose final occurrence the isoform-specific
  # identifier occurs. Sequences sharing the same pre-delimiter sequence are
  # considered isoforms of one another.
  delimiter = sys.argv[2]

  with open(input_fasta_filename) as in_file:
    seq_lengths = calculate_sequence_lengths(in_file, delimiter)
    in_file.seek(0)
    print_longest_isoforms_for_each_gene(in_file, delimiter, seq_lengths)

if __name__ == '__main__':
  main()
