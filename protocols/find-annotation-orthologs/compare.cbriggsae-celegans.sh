#!/bin/bash

function set_env {
  BASEDIR=~/work/jubilant-peanut
  PROTDIR=$BASEDIR/protocols/find-annotation-orthologs
  RUNNAMES=(no-outgroup.blosum{62,80})
}

function create_dirs {
  rm -rf   $BASEDIR/runs/inparanoid.caenorhabditis
  mkdir -p $BASEDIR/runs/inparanoid.caenorhabditis
  cd       $BASEDIR/runs/inparanoid.caenorhabditis

  mkdir proteins

  for foo in ${RUNNAMES[@]}; do
    mkdir $foo
    cp -a ~/.apps/inparanoid/* $foo/
  done

  # Use multiple CPUs.
  for foo in ${RUNNAMES[@]}; do
    sed -i 's/^$blastall = "blastall"/$blastall = "blastall -a28"/' $foo/inparanoid.pl
  done
  for foo in no-outgroup.blosum80; do
    sed -i 's/^$matrix = "BLOSUM62"/$matrix = "BLOSUM80"/' $foo/inparanoid.pl
  done

  cd ../..
}

function filter_isoforms_for_genome {
  id=$1
  organism=$2
  delimiter=$3

  $PROTDIR/filter_isoforms.py $BASEDIR/data/$id/$organism.$id.WS240.protein.fa $delimiter \
    > proteins/$id.filtered_isoforms.fa
}

# Retain only longest sequence for each gene, removing shorter isoforms.
function filter_isoforms {
  cd $BASEDIR/runs/inparanoid.caenorhabditis

  filter_isoforms_for_genome PRJNA10731 c_briggsae '.' &
  filter_isoforms_for_genome PRJNA13758 c_elegans  '-' &

  wait
  cd ../..
}

function munge_fasta {
  protein_path=$1
  id=$2

  cat $protein_path | $PROTDIR/munge_fasta.py $id proteins/$id.munged.fa proteins/$id.name_map.json
}

# Replace protein names with numerical IDs so InParanoid is happy.
function alter_fasta_ids {
  cd $BASEDIR/runs/inparanoid.caenorhabditis

  id=PRJNA13758
  munge_fasta $BASEDIR/data/$id/c_elegans.$id.WS239.protein.fa  $id &
  id=PRJNA10731
  munge_fasta $BASEDIR/data/$id/c_briggsae.$id.WS240.protein.fa $id &

  wait
  cd ../..
}

function run_inparanoid {
  cd $BASEDIR/runs/inparanoid.caenorhabditis

  for foo in no-outgroup.blosum{62,80}; do
    cd $foo
    # InParanoid sucks such that it requires the files it processes to reside
    # in the same directory as the application.
    cp -a ../proteins/{PRJNA13758,PRJNA10731}.munged.fa .
    PATH=~/.apps/blast-2.2.26/bin/:$PATH ./inparanoid.pl PRJNA13758.munged.fa PRJNA10731.munged.fa &
    cd ..
  done

  wait
  cd ../..
}

function process_results {
  cd $BASEDIR/runs/inparanoid.caenorhabditis

  for foo in ${RUNNAMES[@]}; do
    cd $foo
    if [[ ! -d results ]]; then
      mkdir results
    fi

    PARSED=results/inparanoid_results.json
    $PROTDIR/parse_inparanoid.py Output*.fa > $PARSED
    cat $PARSED | $PROTDIR/summarize_groups.py > results/orthologue-summary
    cat $PARSED | $PROTDIR/plot_group_summary.py 'MHco3(ISE)' 'McMaster' results/ortho_groups_{linear,log}.png
    cat $PARSED | $PROTDIR/examine_transcript_groups.py                        \
      ../proteins/{PRJNA13758,PRJNA10731}.name_map.json                         \
      $BASEDIR/data/PRJNA13758/c_elegans.PRJNA13758.WS239.annotations.gff3       \
      $BASEDIR/data/PRJNA10731/c_briggsae.PRJNA10731.WS240.annotations.gff3 \
      > results/transcript-groups

    cd ..
  done

  cd ../..
}

function main {
  set_env
  create_dirs
  alter_fasta_ids
  run_inparanoid
  process_results
}

main
