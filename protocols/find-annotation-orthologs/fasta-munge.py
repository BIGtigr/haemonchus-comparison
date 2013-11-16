#!/usr/bin/env python3
import sys
from dregs import binf

def main():
  count = 1
  for seq_id, seq in binf.parse_fasta(sys.stdin):
    binf.write_fasta_seq(sys.stdout, '%s_prot%s' % (sys.argv[1], count), seq)
    count += 1

main()
