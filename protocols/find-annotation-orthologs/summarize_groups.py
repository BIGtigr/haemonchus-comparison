#!/usr/bin/env python2
'''
Print summary of relationship cardinality in orhologous groups.

Usage: cat inparanoid_results.json | summarize_groups.py
'''

import json
import sys
from collections import defaultdict

def summarize_groups(groups):
  '''
  Suppose InParanoid was run on sequence sets A and B. This function returns a
  dictionary counting occurrences of the following orthologue groupp
  relationships:

    * 1_to_1: Instances where a single sequence in A corresponds to a single
      sequence in B
    * 1_to_n: Instances where a single sequence in A corresponds to multiple
      sequences in B
    * n_to_1: Instances where multiple sequences in A correspond to a single
      sequence in B
    * n_to_n: Instances where multiple sequences in A correspond to multiple
      sequences in B
  '''
  summary = defaultdict(int)

  for group in groups:
    a_len = len(group['a'])
    b_len = len(group['b'])

    if a_len == b_len == 1:
      summary['1_to_1'] += 1
    elif a_len > 1 and b_len > 1:
      summary['n_to_n'] += 1
    elif a_len > 1:
      summary['n_to_1'] += 1
    elif b_len > 1:
      summary['1_to_n'] += 1
    else:
      raise Exception('Unexpectedly empty group')

  return dict(summary)

def main():
  ortho_groups = json.load(sys.stdin)['groups']
  summary = summarize_groups(ortho_groups)

  for key in sorted(summary.keys()):
    print('%s=%s' % (key, summary[key]))

if __name__ == '__main__':
  main()
