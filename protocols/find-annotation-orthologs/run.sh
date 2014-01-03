#!/bin/bash

function set_env {
  BASEDIR=~/work/jubilant-peanut
  PROTDIR=$BASEDIR/protocols/find-annotation-orthologs
  RUNNAMES=({celegans,no}-outgroup.blosum{62,80})
}

function create_dirs {
  rm -rf   $BASEDIR/runs/inparanoid
  mkdir -p $BASEDIR/runs/inparanoid
  cd       $BASEDIR/runs/inparanoid

  for foo in ${RUNNAMES[@]}; do
    mkdir $foo
    cp -a ~/.apps/inparanoid/* $foo/
  done

  # Use multiple CPUs.
  for foo in ${RUNNAMES[@]}; do
    sed -i 's/^$blastall = "blastall"/$blastall = "blastall -a28"/' $foo/inparanoid.pl
  done
  for foo in celegans-outgroup.blosum{62,80}; do
    sed -i 's/^$use_outgroup = 0/$use_outgroup = 1/' $foo/inparanoid.pl
  done
  for foo in {no,celegans}-outgroup.blosum80; do
    sed -i 's/^$matrix = "BLOSUM62"/$matrix = "BLOSUM80"/' $foo/inparanoid.pl
  done

  cd ../..
}

function munge_fasta {
  organism=$1
  id=$2

  cat $BASEDIR/data/$id/$organism.$id.WS239.protein.fa | $PROTDIR/munge_fasta.py $id proteins/$id.munged.fa proteins/$id.name_map.json
}

# Replace protein names with numerical IDs so InParanoid is happy.
function alter_fasta_ids {
  cd $BASEDIR/runs/inparanoid
  mkdir -p proteins/

  munge_fasta h_contortus PRJEB506
  munge_fasta h_contortus PRJNA205202
  munge_fasta c_elegans   PRJNA13758

  cd ../..
}

function run_inparanoid {
  cd $BASEDIR/runs/inparanoid

  for foo in no-outgroup.blosum{62,80}; do
    cd $foo
    # InParanoid sucks such that it requires the files it processes to reside
    # in the same directory as the application.
    cp -a ../proteins/{PRJEB506,PRJNA205202}.munged.fa .
    PATH=~/.apps/blast-2.2.26/bin/:$PATH ./inparanoid.pl PRJEB506.munged.fa PRJNA205202.munged.fa
    cd ..
  done
  for foo in celegans-outgroup.blosum{62,80}; do
    cd $foo
    cp -a ../proteins/{PRJEB506,PRJNA205202,PRJNA13758}.munged.fa .
    PATH=~/.apps/blast-2.2.26/bin/:$PATH ./inparanoid.pl PRJEB506.munged.fa PRJNA205202.munged.fa PRJNA13758.munged.fa
    cd ..
  done

  cd ../..
}

function process_results {
  cd $BASEDIR/runs/inparanoid

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
      ../proteins/{PRJEB506,PRJNA205202}.name_map.json                         \
      $BASEDIR/data/PRJEB506/h_contortus.PRJEB506.WS239.annotations.gff3       \
      $BASEDIR/data/PRJNA205202/h_contortus.PRJNA205202.WS239.annotations.gff3 \
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
