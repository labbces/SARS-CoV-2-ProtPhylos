#!/usr/bin/bash 


#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -t 1-35467




module load   Hmmer/3.2.1 
INFILE=`ls Europe/*.fasta| head -n $SGE_TASK_ID|tail -n 1`
BASENAME=${INFILE/\.fasta}
GFF=${BASENAME}.gff
prodigal  -i $INFILE -o ${BASENAME}.gff  -f gff -a ${BASENAME}.proteins.faa -d ${BASENAME}.nucleotides.fna -p meta



