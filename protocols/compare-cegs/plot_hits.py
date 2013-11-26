import json
import sys
from pprint import pprint
from collections import Counter
from functools import reduce

def parse_hits(hits_filename):
  with open('hits.json') as hits_file:
    raw_hits = hits_file.read()
  hits = json.loads(raw_hits)
  return hits

def calc_counts(sets):
  if len(sets) != 4:
    raise Exception('Need four sets')
  sets = [set(s) for s in sets]
  union = reduce(lambda x, y: x.union(y), sets)

  indicators = [4*'%d' % (a in sets[0], a in sets[1], a in sets[2], a in sets[3]) for a in union]
  counts = Counter(indicators)
  all_counts = {}
  for a in (0, 1):
    for b in (0, 1):
      for c in (0, 1):
        for d in (0, 1):
          key = 4*'%s' % (a, b, c, d)
          all_counts[key] = counts[key]
  return all_counts

def plot(labels, counts):
  labels_str = ', '.join(["'%s'" % s for s in labels])
  weights_str = ', '.join(["'%s' = %s" % (k, v) for k, v in counts.items()])

  print('''
    library(Vennerable)
    v <- Venn(SetNames = c(%s), Weight = c(%s))
    svg('mhco3_vs_mcmaster_cegs.svg')
    plot(v, doWeights = FALSE, type='ellipses', show = list(SetLabels=TRUE))
    dev.off()
  ''' % (labels_str, weights_str))

def main():
  hits = parse_hits(sys.argv[1])
  sources = ('PRJEB506', 'PRJNA205202')

  labels = ('MHco3 complete', 'MHco3 partial', 'McMaster complete', 'McMaster partial')
  sets = (
    hits[sources[0]]['complete'],
    hits[sources[0]]['partial'],
    hits[sources[1]]['complete'],
    hits[sources[1]]['partial'],
  )

  counts = calc_counts(sets)
  plot(labels, counts)

main()
