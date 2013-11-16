import json
import os

def main():
  hits_json = open('hits.json').read()
  hits = json.loads(hits_json)

  for i in hits.keys():
    for j in hits[i].keys():
      hits[i][j] = set(hits[i][j])

  source_sets = ('PRJEB506', 'PRJNA205202')
  write_missing(hits, *source_sets)
  write_missing(hits, *reversed(source_sets))
  write_complete_vs_partial(hits, *source_sets)
  write_complete_vs_partial(hits, *reversed(source_sets))

def write_set(happy_set, fname):
  with open(os.path.join('ceg-comparisons', fname), 'w') as fout:
    fout.write('\n'.join(happy_set) + '\n')

def write_missing(hits, name_a, name_b):
  missing = hits[name_a]['complete'].union(hits[name_a]['partial']) - hits[name_b]['complete'].union(hits[name_b]['partial'])
  write_set(missing, '%s_minus_%s.set' % (name_a, name_b))

def write_complete_vs_partial(hits, name_a, name_b):
  intersection = hits[name_a]['complete'].intersection(hits[name_b]['partial'])
  write_set(intersection, '%s_complete_intersection_%s_partial.set' % (name_a, name_b))

main()
