#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 2
#$ -t 1-3

FASTA=`ls -1 nsp*.fasta | head -n $SGE_TASK_ID | tail -n  1`
HMM=${FASTA/\.fasta/_single.hmm}
module load Hmmer/3.2.1
hmmbuild --amino --cpu $NSLOTS --seed 26 $HMM $FASTA
