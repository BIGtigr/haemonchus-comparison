#!/usr/bin/env python2

import json
import re
import sqlite3
import sys

class TranscriptManager(object):
  def __init__(self):
    self._conn = sqlite3.connect(':memory:')
    # Give dictionary-based access to results.
    self._conn.row_factory = sqlite3.Row
    self._create_db()

  def _create_db(self):
    cursor = self._cursor()
    cursor.execute('''CREATE TABLE transcripts (
      id         INTEGER PRIMARY KEY,
      group_name TEXT,
      name       TEXT,
      seqid      TEXT,
      start      INTEGER,
      end        INTEGER,
      strand     TEXT
    )''')
    cursor.execute('CREATE INDEX idx_name ON transcripts(name)')
    cursor.execute('CREATE INDEX idx_seqid_strand ON transcripts(seqid, strand)')
    self._commit()

  def _insert_transcripts(self, transcripts, group_name):
    cursor = self._cursor()

    for tid, transcript in transcripts.items():
      cursor.execute(
        'INSERT INTO transcripts (group_name, name, seqid, start, end, strand) VALUES (?, ?, ?, ?, ?, ?)', (
        group_name,
        tid,
        transcript['seqid'],
        transcript['start'],
        transcript['end'],
        transcript['strand'],
      ))

    self._commit()

  def _parse_gff(self, gff_filename):
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
        tid = attribs['ID']
        if tid in transcripts:
          raise Exception('Duplicate seqid %s' % seqid)
        transcripts[tid] = {
          'seqid':  fields[0],
          'start':  int(fields[3]),
          'end':    int(fields[4]),
          'strand': fields[6],
        }

    return transcripts

  def _cursor(self):
    return self._conn.cursor()

  def _commit(self):
    self._conn.commit()

  def parse(self, gff_filename, group_name):
    transcripts = self._parse_gff(gff_filename)
    self._insert_transcripts(transcripts, group_name)

  def close(self):
    self._conn.close()

  def find_scaffolds(self, transcript_names):
    cursor = self._cursor()
    placeholders = ', '.join(['?' for t in transcript_names])
    cursor.execute('''
      SELECT
        t1.group_name,
        t1.seqid,
        t1.strand,
        COUNT(t1.seqid) AS orthologues_on_scaffold_count,
        (
          SELECT COUNT(*)
          FROM transcripts t2
          WHERE
            t2.seqid = t1.seqid AND
            t2.strand = t1.strand
        ) AS transcripts_on_scaffold_count
      FROM transcripts t1
      WHERE t1.name in (%s)
      GROUP BY t1.seqid, t1.strand, t1.group_name
      ORDER BY t1.group_name, t1.seqid
      ''' % placeholders,
      ['transcript:' + tname for tname in transcript_names]
    )
    rows = cursor.fetchall()

    # As any transcripts not in DB for "WHERE t1.name IN (...)" will silently
    # be omitted from SQL result set, we must ensure that the sizes of the
    # query and results sets match.
    transcripts_on_scaffolds_sum = sum([t['orthologues_on_scaffold_count'] for t in rows])
    assert transcripts_on_scaffolds_sum == len(transcript_names), \
      'Some transcripts lack associated scaffolds'

    return rows


# Natural sorting technique taken from http://stackoverflow.com/a/2669120/1691611.
#def sort_scaffold_counts(scaffold_counts):
  #convert      = lambda text: int(text) if text.isdigit() else text
  #alphanum_key = lambda key:  [convert(c) for c in re.split('([0-9]+)', key)]
  #scaffold_counts.sort(key=lambda c: alphanum_key(c[1]))

def process_ortho_group(ogroup, mappings, transcript_manager):
  for gname in ogroup.keys():
    transcripts = []
    for seqname, score in ogroup[gname]:
      orig_full_name = mappings[gname][seqname]
      orig_name      = orig_full_name.split()[0]
      transcripts.append(orig_name)

    for row in transcript_manager.find_scaffolds(transcripts):
      print('%s %s %s %s %s' % (
        row['group_name'],
        row['seqid'].ljust(20),
        row['strand'].ljust(4),
        row['orthologues_on_scaffold_count'],
        row['transcripts_on_scaffold_count'],
      ))

def process_ortho_groups(ortho_groups, mappings, transcript_manager):
  already_printed_first = False

  for i, ogroup in enumerate(ortho_groups):
    a_len, b_len = len(ogroup['a']), len(ogroup['b'])
    if a_len == b_len == 1:
      continue

    if already_printed_first:
      print('')
    else:
      already_printed_first = True

    print('Group (a:b = %s:%s)' % (a_len, b_len))
    process_ortho_group(ogroup, mappings, transcript_manager)

def examine_contiguity(ortho_groups, transcript_mapping_fnames, transcript_fnames):
  tm       = TranscriptManager()
  mappings = {}

  for gname in ('a', 'b'):
    with open(transcript_mapping_fnames[gname]) as f:
      mappings[gname] = json.load(f)
    tm.parse(transcript_fnames[gname], gname)

  process_ortho_groups(ortho_groups, mappings, tm)

def main():
  transcript_mapping_fnames = {'a': sys.argv[1], 'b': sys.argv[2]}
  transcript_fnames         = {'a': sys.argv[3], 'b': sys.argv[4]}

  ortho_groups = json.load(sys.stdin)['groups']
  examine_contiguity(ortho_groups, transcript_mapping_fnames, transcript_fnames)

if __name__ == '__main__':
  #import cProfile
  #cProfile.run('main()')
  main()
