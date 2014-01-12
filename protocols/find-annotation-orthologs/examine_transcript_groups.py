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
      tstart     INTEGER,
      tend       INTEGER,
      strand     TEXT
    )''')
    cursor.execute('CREATE INDEX idx_name ON transcripts(name)')
    cursor.execute('CREATE INDEX idx_seqid_strand ON transcripts(seqid, strand)')
    self._commit()

  def _insert_transcripts(self, transcripts, group_name):
    cursor = self._cursor()

    for tid, transcript in transcripts.items():
      cursor.execute(
        'INSERT INTO transcripts (group_name, name, seqid, tstart, tend, strand) VALUES (?, ?, ?, ?, ?, ?)', (
        group_name,
        tid,
        transcript['seqid'],
        transcript['tstart'],
        transcript['tend'],
        transcript['strand'],
      ))

    self._commit()

  def _retrieve_transcript_names(self, mapping):
    '''
    Retrieve transcript names specified in mapping file. This is necessary as
    not all transcripts listed in the GFF file should be inserted into the
    database -- we want only the transcripts specified in the input FASTA file.
    (Note that, in the current implementation, this FASTA file will list only
    the longest isoform for each gene.)
    '''
    gene_names = mapping.values()
    gene_names = [n.split()[0] for n in gene_names]
    return gene_names

  def _parse_gff(self, gff_filename, mapping):
    transcripts = {}
    mapped_transcript_names = self._retrieve_transcript_names(mapping)

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
        if tid.startswith('transcript:'):
          tid = tid.split(':', 1)[1]

        if tid in transcripts:
          raise Exception('Duplicate seqid %s' % seqid)
        if tid not in mapped_transcript_names:
          continue

        transcripts[tid] = {
          'seqid':  fields[0],
          'tstart':  int(fields[3]),
          'tend':    int(fields[4]),
          'strand': fields[6],
        }

    return transcripts

  def _cursor(self):
    return self._conn.cursor()

  def _commit(self):
    self._conn.commit()

  def parse(self, gff_filename, mapping, group_name):
    transcripts = self._parse_gff(gff_filename, mapping)
    self._insert_transcripts(transcripts, group_name)

  def close(self):
    self._conn.close()

  def _determine_intervening_genes(self, row, placeholders, transcript_names):
    cursor = self._cursor()
    cursor.execute('''
      SELECT COUNT(*) AS intervening_genes
      FROM transcripts t1
      WHERE
        t1.name NOT IN (%s)  AND
        t1.group_name = ?    AND
        t1.seqid = ?         AND
        t1.strand = ?        AND
        (
          (CASE WHEN t1.tend   < ? THEN t1.tend   ELSE ? END) -- MIN(t1.tend, placeholder)
          >
          (CASE WHEN t1.tstart > ? THEN t1.tstart ELSE ? END) -- MAX(t1.tstart, placeholder)
        )
      ''' % placeholders, transcript_names + [
        row['group_name'],
        row['seqid'],
        row['strand'],
        row['group_end'],
        row['group_end'],
        row['group_start'],
        row['group_start'],
      ]
    )

    result = cursor.fetchone()
    row['intervening_genes'] = result['intervening_genes']

  def find_scaffolds(self, transcript_names):
    cursor = self._cursor()
    placeholders     = ', '.join(['?' for t in transcript_names])
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
        ) AS transcripts_on_scaffold_count,
        MIN(t1.tstart) AS group_start,
        MAX(t1.tend)   AS group_end
      FROM transcripts t1
      WHERE t1.name IN (%s)
      GROUP BY t1.seqid, t1.strand, t1.group_name
      ORDER BY t1.group_name, t1.seqid
      ''' % placeholders,
      transcript_names
    )
    rows = cursor.fetchall()
    rows = [dict(r) for r in rows]
    for row in rows:
      self._determine_intervening_genes(row, placeholders, transcript_names)

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
      print('%s %-20s %-4s %2s %2s %-5s %s' % (
        row['group_name'],
        row['seqid'],
        row['strand'],
        row['orthologues_on_scaffold_count'],
        row['transcripts_on_scaffold_count'],
        # Boolean indicating whether orthologues tandemly arrayed.
        row['intervening_genes'] == 0,
        row['intervening_genes'],
      ))

def process_ortho_groups(ortho_groups, mappings, transcript_manager):
  already_printed_first = False

  for i, ogroup in enumerate(ortho_groups):
    a_len, b_len = len(ogroup['a']), len(ogroup['b'])
    # We aren't interested in 1:1 relationships between orthologues, as we
    # already know each orthologue lies on a single scaffold.
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
    tm.parse(transcript_fnames[gname], mappings[gname], gname)

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
