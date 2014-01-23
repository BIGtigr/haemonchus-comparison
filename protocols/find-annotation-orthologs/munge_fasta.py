#!/usr/bin/env python3

import json
import sys
from dregs import binf

def main():
  seq_set_id = sys.argv[1]
  munged_fasta_filename = sys.argv[2]
  mapping_filename = sys.argv[3]

  name_mapping = {}
  count = 1

  with open(munged_fasta_filename, 'w') as munged_fasta_file:
    for seq_id, seq in binf.parse_fasta(sys.stdin):
      new_name = '%s_prot%s' % (seq_set_id, count)
      name_mapping[new_name] = seq_id
      binf.write_fasta_seq(munged_fasta_file, new_name, seq)
      count += 1

  with open(mapping_filename, 'w') as mapping_file:
    json.dump(name_mapping, mapping_file)

main()
