#!/bin/bash
OUTDIR=../../runs/inparanoid
mkdir -p $OUTDIR
#EXTRA_ARGS="--only-analyze"

for matrix in BLOSUM{62,80}; do
  # Haemonchus against Haemonchus.
  ./run_inparanoid.py \
    --matrix $matrix \
    --outgroup PRJNA13758 \
    $EXTRA_ARGS \
    PRJEB506 PRJNA205202 \
    $OUTDIR &
  ./run_inparanoid.py \
    --matrix $matrix \
    $EXTRA_ARGS \
    PRJEB506 PRJNA205202 \
    $OUTDIR &

  # C. elegans against C. briggsae.
  ./run_inparanoid.py \
    --matrix $matrix \
    $EXTRA_ARGS \
    PRJNA13758 PRJNA10731 \
    $OUTDIR &

  # Haemonchus against C. elegans.
  ./run_inparanoid.py \
    --matrix $matrix \
    $EXTRA_ARGS \
    PRJEB506 PRJNA13758 \
    $OUTDIR &
  ./run_inparanoid.py \
    --matrix $matrix \
    $EXTRA_ARGS \
    PRJNA205202 PRJNA13758 \
    $OUTDIR &
done

wait
