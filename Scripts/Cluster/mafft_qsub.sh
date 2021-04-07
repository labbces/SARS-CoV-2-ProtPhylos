#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -pe smp 20

module load mafft/7.407  

mafft --thread $NSLOTS --auto Africa_Protein+E_sem_redu.fasta > Africa_Protein+E_sem_redu.aln

   



