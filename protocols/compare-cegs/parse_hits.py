import json
import sys
from collections import defaultdict
from pprint import pprint

def parse_hits(hits_filename):
  hits = defaultdict(lambda: defaultdict(list))

  with open(hits_filename) as hits_file:
    for line in hits_file:
      tokens = line.strip().split()
      fields = dict([t.split('=', 1) for t in tokens])

      hit_type = fields['type'].lower()
      kog      = fields['kog']
      del fields['type'], fields['kog']

      for k, v in fields.items():
        try:
          fields[k] = int(v)
        except ValueError:
          try:
            fields[k] = float(v)
          except ValueError:
            pass

      hits[hit_type][kog].append(fields)

  return hits

def extract_kog_names(hits):
  names = {}
  for hit_type in hits.keys():
    names[hit_type] = set(hits[hit_type].keys())
    # Convert to list to permit JSON serialization
    names[hit_type] = list(names[hit_type])
  return names

def find_overlap(hit_names):
  for hit_type, hits in hit_names.items():
    pass

def remove_full_hits_from_partials(hits):
  hit_types = set(hits.keys())
  if set(('complete', 'partial')) != hit_types:
    raise Exception('Unexpected hit types: %s' % hit_types)

  for kog in hits['complete'].keys():
    if kog in hits['partial'].keys():
      del hits['partial'][kog]

def main():
  all_hits = {}

  for hit_source in sys.argv[1:]:
    hit_name, hit_filename = hit_source.split('=', 1)

    hits_for_file = parse_hits(hit_filename)
    remove_full_hits_from_partials(hits_for_file)
    hit_names = extract_kog_names(hits_for_file)
    all_hits[hit_name] = hit_names
  print(json.dumps(all_hits))

main()
