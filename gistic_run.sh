#!/bin/bash

#### GISTIC2 copy number significance for Canine (Dog) CanFam 3.1 genome
## https://github.com/sbamin/canine_gistic2

## Authors:
# Samir B. Amin (@sbamin)
# Emmanuel Martinez (@jemartinezledes)

## load MCR/8.0 env; other versions may not work
module load rvMCR/8.0
sleep 2
echo "MCR_ROOT is ${MCR_ROOT}"

## set vars
## gistic code directory
GITDIR="/projects/verhaak-lab/verhaak_env/verhaak_apps/canine_gistic2"
## output dir
BASEDIR="/scratch/output"
## canfam3_1_order.mat was built by @jemartinezledes
REFGENE="${GITDIR}/canfam3_1_order.mat"

# export GISTIC2 code directory in PATH
export PATH="${PATH}:${GITDIR}"

## run with scna segment file
OUTDIR1="${BASEDIR}/case1_viterbi"
mkdir -p "${OUTDIR1}"
## Read README on how to make segment and marker file.
SEGFILE1="${BASEDIR}/inputs/gistic2_segments.tsv"
MARKERS1="${BASEDIR}/inputs/markers_gistic.txt"

#### RUN GISTIC2 ####
cd "$GITDIR" && \
gp_gistic2_from_seg_upd -b "${OUTDIR1}" -seg "${SEGFILE1}" -mk "${MARKERS1}" -refgene "${REFGENE}" -maxseg 25000 -savegene 1 -genegistic 1
echo "GISTIC2 run ended with exit code $?"

## END ##
