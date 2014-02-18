#!/usr/bin/env python3

import argparse
import glob
import sarge
import os

class CommandManager(object):
  def __init__(self):
    self._pipelines = []

  def run(self, *args, **kwargs):
    print(args)
    p = sarge.run(*args, **kwargs)
    if 'async' in kwargs and kwargs['async']:
      self._pipelines.append(p)

  def wait(self):
    for p in self._pipelines:
      p.close()
    self._pipelines = []

CONF = {}
cm = CommandManager()

def set_config(run_dir):
  global CONF
  CONF['BASEDIR'] = os.path.abspath(os.path.expanduser('~/work/jubilant-peanut'))
  CONF['PROTDIR'] = CONF['BASEDIR'] + '/protocols/find-annotation-orthologs'
  CONF['RUNDIR']  = os.path.abspath(run_dir)

  CONF['DATASETS'] = {
    'PRJEB506': {
      'proteins':   CONF['BASEDIR'] + '/data/PRJEB506/h_contortus.PRJEB506.WS239.protein.fa',
      'annotation': CONF['BASEDIR'] + '/data/PRJEB506/h_contortus.PRJEB506.WS239.annotations.gff3',
      'filter_isoforms': True,
      'isoform_delimiter': '.',
      'friendly_name': 'MHco3(ISE)',
    },
    'PRJNA205202': {
      'proteins':   CONF['BASEDIR'] + '/data/PRJNA205202/h_contortus.PRJNA205202.WS239.protein.fa',
      'annotation': CONF['BASEDIR'] + '/data/PRJNA205202/h_contortus.PRJNA205202.WS239.annotations.gff3',
      'filter_isoforms': True,
      'isoform_delimiter': '-',
      'friendly_name': 'McMaster',
    },
    'PRJNA13758': {
      'proteins':   CONF['BASEDIR'] + '/data/PRJNA13758/c_elegans.PRJNA13758.WS239.protein.fa',
      'annotation': CONF['BASEDIR'] + '/data/PRJNA13758/c_elegans.PRJNA13758.WS239.annotations.gff3',
      'filter_isoforms': False,
      'friendly_name': 'C. elegans',
    },
    'PRJNA10731': {
      'proteins':   CONF['BASEDIR'] + '/data/PRJNA10731/c_briggsae.PRJNA10731.WS240.protein.fa',
      'annotation': CONF['BASEDIR'] + '/data/PRJNA10731/c_briggsae.PRJNA10731.WS240.annotations.gff3',
      'filter_isoforms': False,
      'friendly_name': 'C. briggsae',
    },
  }

def create_dirs():
  cm.run('rm -rf ' +  CONF['RUNDIR'])
  cm.run('mkdir -p ' + CONF['RUNDIR'])
  os.chdir(CONF['RUNDIR'])
  cm.run('mkdir proteins')

def filter_isoforms(genome_id, protein_path, delimiter):
  filtered_file = '{rundir}/proteins/{genome_id}.filtered_isoforms.fa'.format(
    rundir = CONF['RUNDIR'],
    genome_id = genome_id,
  )
  cm.run('{protdir}/filter_isoforms.py {protein_path} {delimiter} > {filtered_file}'.format(
    protdir = CONF['PROTDIR'],
    genome_id = genome_id,
    protein_path = protein_path,
    delimiter = delimiter,
    filtered_file = filtered_file,
  ))
  return filtered_file

def alter_fasta_ids(genome_id, protein_path):
  munged_file = '{rundir}/proteins/{genome_id}.munged.fa'.format(
    rundir = CONF['RUNDIR'],
    genome_id = genome_id,
  )
  name_map = '{rundir}/proteins/{genome_id}.name_map.json'.format(
    rundir = CONF['RUNDIR'],
    genome_id = genome_id,
  )

  cm.run('cat {protein_path} | {protdir}/munge_fasta.py {genome_id} {munged_file} {name_map}'.format(
    protein_path = protein_path,
    protdir = CONF['PROTDIR'],
    genome_id = genome_id,
    munged_file = munged_file,
    name_map = name_map,
  ))
  return {
    'munged_fasta': munged_file,
    'name_map':     name_map
  }


def prepare_inparanoid(run_name, outgroup):
  os.chdir(CONF['RUNDIR'])
  cm.run('cp -a {inparanoid_path} {run_name}'.format(
    run_name = run_name,
    inparanoid_path = os.path.expanduser('~/.apps/inparanoid/'),
  ))

  if outgroup:
    cm.run('''sed -i 's/^$use_outgroup = 0/$use_outgroup = 1/' {run_name}/inparanoid.pl'''.format(run_name = run_name))
  # Use multiple CPUs.
  cm.run('''sed -i 's/^$blastall = "blastall"/$blastall = "blastall -a28"/' {run_name}/inparanoid.pl'''.format(run_name = run_name))
  # Set matrix.
  cm.run('''sed -i 's/^$matrix = "BLOSUM62"/$matrix = "BLOSUM80"/' {run_name}/inparanoid.pl'''.format(run_name = run_name))

def run_inparanoid(run_name, genome_a, genome_b, outgroup):
  os.chdir(os.path.join(CONF['RUNDIR'], run_name))

  genomes = [genome_a, genome_b]
  if outgroup:
    genomes.append(outgroup)
  fasta_files = ['%s.munged.fa' % g for g in genomes]

  for fasta_file in fasta_files:
    cm.run('cp -a ../proteins/{fasta_file} .'.format(fasta_file = fasta_file))

  cm.run(
    './inparanoid.pl {fasta_files}'.format(fasta_files = ' '.join(fasta_files)),
    env = {'PATH': '{blast_path}:{existing_path}'.format(
      blast_path = os.path.expanduser('~/.apps/blast-2.2.26/bin'),
      existing_path = os.environ['PATH'],
    )}
  )

def process_results(run_name, genome_a, genome_b, munged_fnames):
  os.chdir(os.path.join(CONF['RUNDIR'], run_name))
  cm.run('mkdir results')
  parsed = 'results/inparanoid_results.json'

  cm.run('{protdir}/parse_inparanoid.py {inparanoid_output} > {parsed}'.format(
    protdir = CONF['PROTDIR'],
    inparanoid_output = glob.glob('Output*')[0],
    parsed = parsed
  ))
  cm.run('cat {parsed} | {protdir}/summarize_groups.py > results/orthologue-summary'.format(
    parsed = parsed,
    protdir = CONF['PROTDIR']
  ))
  cm.run(("cat {parsed} | {protdir}/plot_group_summary.py '{friendly_name_a}' '{friendly_name_b}' " + \
    "results/ortho_groups_linear.png results/ortho_groups_log.png").format(
    parsed = parsed,
    protdir = CONF['PROTDIR'],
    friendly_name_a = CONF['DATASETS'][genome_a]['friendly_name'],
    friendly_name_b = CONF['DATASETS'][genome_b]['friendly_name'],
  ))
  cm.run(('cat {parsed} | {protdir}/examine_transcript_groups.py ' + \
    '{name_map_a} ' + \
    '{name_map_b} ' + \
    '{annotation_a} ' + \
    '{annotation_b} ' + \
    '> results/transcript-groups').format(
    parsed = parsed,
    protdir = CONF['PROTDIR'],
    name_map_a = munged_fnames[genome_a]['name_map'],
    name_map_b = munged_fnames[genome_b]['name_map'],
    annotation_a = CONF['DATASETS'][genome_a]['annotation'],
    annotation_b = CONF['DATASETS'][genome_b]['annotation'],
  ))

def perform_run(genome_a, genome_b, outgroup, matrix, munged_fnames):
  run_name = '{a}-{b}.{outgroup}-outgroup.{matrix}'.format(
    a = genome_a,
    b = genome_b,
    outgroup = outgroup or 'no',
    matrix = matrix,
  )
  prepare_inparanoid(run_name, outgroup)
  run_inparanoid(run_name, genome_a, genome_b, outgroup)
  process_results(run_name, genome_a, genome_b, munged_fnames)

# This could be moved into a different function, as we only need to munge the
# FASTA files for genomes used in our current run, rather than all specified
# within CONFIG['DATASETS']. Doing so would only reduce the number of genomes
# to munge from four to three (when an outgroup is used), however, which would
# reduce the runtime by only a few minutes. Thus, I choose instead to save
# myself the effort of restructuring this operation.
def munge_all_fasta_files(datasets):
  munged = {}

  for genome_id, params in datasets.items():
    protein_path = params['proteins']
    if params['filter_isoforms']:
      protein_path = filter_isoforms(genome_id, protein_path, params['isoform_delimiter'])
    munged[genome_id] = alter_fasta_ids(genome_id, params['proteins'])

  return munged

def main():
  parser = argparse.ArgumentParser(description='Run InParanoid and parse results.')
  parser.add_argument('genome_a',   help='First genome')
  parser.add_argument('genome_b',   help='Second genome')
  parser.add_argument('output_dir', help='Output directory')
  parser.add_argument('-o', '--outgroup', dest='outgroup', action='store', help='Outgroup used to break up large orthologous groups')
  parser.add_argument('-m', '--matrix', dest='matrix', action='store',
    choices=('BLOSUM62', 'BLOSUM80'), default='BLOSUM62', help='Scoring matrix')
  args = parser.parse_args()

  set_config(args.output_dir)
  create_dirs()
  munged_fnames = munge_all_fasta_files(CONF['DATASETS'])

  perform_run(args.genome_a, args.genome_b, args.outgroup, args.matrix, munged_fnames)

if __name__ == '__main__':
  main()
