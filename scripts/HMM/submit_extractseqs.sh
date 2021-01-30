#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -t 1-4
#$ -pe smp 5

module load EMBOSS/6.6.0 
module load CD-HIT/4.8.1
module load mafft/7.407  
module load Hmmer/3.2.1
  
PROTS=('NSP2' 'NSP6' 'NSP11')
ORDER=$((SGE_TASK_ID-1))
MY_PROT=${PROTS[$ORDER]}
seqret allprot1021.fixnames.fasta:${MY_PROT}_* ${MY_PROT}_all.fasta

cd-hit -i ${MY_PROT}_all.fasta -o ${MY_PROT}_allnr100.fasta -c 1 -M 0 -T $NSLOTS -d 0  

mafft --thread $NSLOTS --auto ${MY_PROT}_allnr100.fasta > ${MY_PROT}_allnr100.aln

hmmbuild --amino --cpu $NSLOTS --seed 26 ${MY_PROT}_allnr100.hmm ${MY_PROT}_allnr100.aln

