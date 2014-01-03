#!/usr/bin/env python2

import json
import re
import sys

def parse_transcripts(gff_filename):
  transcripts = {}

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

# Natural sorting technique taken from http://stackoverflow.com/a/2669120/1691611.
def sort_transcript_tuples(transcript_tuples):
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
  transcript_tuples.sort(key=lambda c: alphanum_key(c[1]))

def examine_contiguity(groups, a_transcript_mappings, b_transcript_mappings, a_transcripts, b_transcripts):
  mappings = {}
  with open(a_transcript_mappings) as f:
    mappings['a'] = json.load(f)
  with open(b_transcript_mappings) as f:
    mappings['b'] = json.load(f)

  transcripts = {
    'a': parse_transcripts(a_transcripts),
    'b': parse_transcripts(b_transcripts),
  }

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

    for gname in group.keys():
      group_seqs = []
      for seqname, score in group[gname]:
        seq = {
          'full_name': mappings[gname][seqname],
          'renamed':   seqname,
          'score':     score,
        }
        seq['name'] = seq['full_name'].split()[0]
        seq['transcript'] = transcripts[gname]['transcript:%s' % seq['name']]
        group_seqs.append(seq)

      transcript_names = [(t['transcript']['seqid'], t['transcript']['strand']) for t in group_seqs]
      counts = [(transcript_names.count(t), t[0], t[1]) for t in set(transcript_names)]
      sort_transcript_tuples(counts)

      for count, seqid, strand in counts:
        print('%s %s %s %s' % (gname, seqid.ljust(20), strand.ljust(4), count))

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
