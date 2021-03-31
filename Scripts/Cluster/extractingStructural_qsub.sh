#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -t 1-2921

module load Python/3.7.2

FASTA=`ls -1 ../*.embl | head -n $SGE_TASK_ID | tail -n 1`

time python3  ../../../../Scripts/Cluster/extractingStructural.py  ${FASTA} 

