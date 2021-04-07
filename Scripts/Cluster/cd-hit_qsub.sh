#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -pe smp 20

module load CD-HIT/4.8.1

cd-hit -i Africa_Protein+S.fasta -o Africa_Protein+S_sem_redu.fasta -c 1 -M 0 -T $NSLOTS -d 0

   



