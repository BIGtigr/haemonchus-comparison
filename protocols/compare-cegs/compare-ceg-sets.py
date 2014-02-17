import json
import os
import sys

def write_set(happy_set, fname):
  with open(fname, 'w') as fout:
    fout.write('\n'.join(happy_set) + '\n')

def write_subtraction(hits, name_a, name_b):
  missing = hits[name_a]['complete'].union(hits[name_a]['partial']) - hits[name_b]['complete'].union(hits[name_b]['partial'])
  write_set(missing, '%s_minus_%s.set' % (name_a, name_b))

def write_complete_vs_partial(hits, name_a, name_b):
  intersection = hits[name_a]['complete'].intersection(hits[name_b]['partial'])
  write_set(intersection, '%s_complete_intersection_%s_partial.set' % (name_a, name_b))

def write_cegs_missing_from_both(hits, name_a, name_b, all_cegs_fname, include_partial):
  all_cegs = open(all_cegs_fname).readlines()
  all_cegs = set([l.strip() for l in all_cegs])

  a_hits = hits[name_a]['complete']
  b_hits = hits[name_b]['complete']
  if include_partial:
    a_hits = a_hits.union(hits[name_a]['partial'])
    b_hits = b_hits.union(hits[name_b]['partial'])

  missing_from_a = all_cegs - a_hits
  missing_from_b = all_cegs - b_hits
  missing_from_both = missing_from_a.intersection(missing_from_b)

  hit_type = include_partial and 'comp+part' or 'complete'
  write_set(missing_from_a, 'missing_%s_from_%s.set' % (hit_type, name_a))
  write_set(missing_from_b, 'missing_%s_from_%s.set' % (hit_type, name_b))
  write_set(missing_from_both, 'missing_%s_from_%s_and_%s.set' % (hit_type, name_a, name_b))

def main():
  with open(sys.argv[1]) as hits_file:
    hits = json.load(hits_file)
  for i in hits.keys():
    for j in hits[i].keys():
      hits[i][j] = set(hits[i][j])

  source_sets = ('PRJEB506', 'PRJNA205202')
  write_subtraction(hits, *source_sets)
  write_subtraction(hits, *reversed(source_sets))
  write_complete_vs_partial(hits, *source_sets)
  write_complete_vs_partial(hits, *reversed(source_sets))
  write_cegs_missing_from_both(hits, source_sets[0], source_sets[1], sys.argv[2], False)
  write_cegs_missing_from_both(hits, source_sets[0], source_sets[1], sys.argv[2], True)

main()
