#!/usr/bin/env python2

import json
import re
import sys
from collections import defaultdict

def parse_transcripts(gff_filename):
  transcripts = {}
  transcript_counts = {}

  with open(gff_filename) as gff_file:
    for line in gff_file:
      line = line.strip()
      if line == '' or line.startswith('#'):
        continue

      fields = line.split('\t')
      if not fields[2] == 'mRNA':
        continue

      attribs = dict([f.split('=', 1) for f in fields[8].split(';')])
      seqid = attribs['ID']
      if seqid in transcripts:
        raise Exception('Duplicate seqid %s' % seqid)
      transcripts[seqid] = {
        'seqid': fields[0],
        'start': int(fields[3]),
        'end':  int(fields[4]),
        'strand': fields[6],
      }

  return transcripts

def count_transcripts_per_seq(transcripts):
  counts = defaultdict(int)
  for transcript in transcripts.values():
    key = (transcript['seqid'], transcript['strand'])
    counts[key] += 1
  return counts

# Natural sorting technique taken from http://stackoverflow.com/a/2669120/1691611.
def sort_scaffold_counts(scaffold_counts):
  convert      = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key:  [convert(c) for c in re.split('([0-9]+)', key)]
  scaffold_counts.sort(key=lambda c: alphanum_key(c[1]))

def process_group(group, mappings, transcripts, transcript_counts_on_scaffolds):
  for gname in group.keys():
    genes = []
    for seqname, score in group[gname]:
      gene = {
        'full_name': mappings[gname][seqname],
        'renamed':   seqname,
        'score':     score,
      }
      gene['name']       = gene['full_name'].split()[0]
      gene['transcript'] = transcripts[gname]['transcript:%s' % gene['name']]
      genes.append(gene)

    transcript_names = [(g['transcript']['seqid'], g['transcript']['strand']) for g in genes]
    scaffold_counts  = [(transcript_names.count(t), t[0], t[1]) for t in set(transcript_names)]
    sort_scaffold_counts(scaffold_counts)

    for count, seqid, strand in scaffold_counts:
      print('%s %s %s %s %s' % (
        gname,
        seqid.ljust(20),
        strand.ljust(4),
        count,
        transcript_counts_on_scaffolds[gname][(seqid, strand)],
      ))

def process_groups(groups, mappings, transcripts, transcript_counts):
  already_printed_first = False

  for i, group in enumerate(groups):
    a_len, b_len = len(group['a']), len(group['b'])
    if a_len == b_len == 1:
      continue

    if already_printed_first:
      print('')
    else:
      already_printed_first = True

    print('Group (a:b = %s:%s)' % (a_len, b_len))
    process_group(group, mappings, transcripts, transcript_counts)

def examine_contiguity(groups, a_transcript_mappings, b_transcript_mappings, a_transcripts, b_transcripts):
  transcript_mapping_filenames = {
    'a': a_transcript_mappings,
    'b': b_transcript_mappings,
  }
  transcript_filenames = {
    'a': a_transcripts,
    'b': b_transcripts,
  }

  mappings          = {}
  transcripts       = {}
  transcript_counts = {}

  for gname in ('a', 'b'):
    with open(transcript_mapping_filenames[gname]) as f:
      mappings[gname]        = json.load(f)
    transcripts[gname]       = parse_transcripts(transcript_filenames[gname])
    transcript_counts[gname] = count_transcripts_per_seq(transcripts[gname])

  process_groups(groups, mappings, transcripts, transcript_counts)

def main():
  group_a_transcript_mappings = sys.argv[1]
  group_b_transcript_mappings = sys.argv[2]
  group_a_transcripts         = sys.argv[3]
  group_b_transcripts         = sys.argv[4]

  ortho_groups = json.load(sys.stdin)['groups']
  examine_contiguity(
    ortho_groups,
    group_a_transcript_mappings,
    group_b_transcript_mappings,
    group_a_transcripts,
    group_b_transcripts,
  )

if __name__ == '__main__':
  main()
