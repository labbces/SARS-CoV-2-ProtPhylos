#!/usr/bin/bash 


#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -t 1-593


module load Prodigal/2.6.3
INFILE=`ls Africa/*.fasta| head -n $SGE_TASK_ID|tail -n 1`
BASENAME=${INFILE/\.fasta}
GFF=${BASENAME}.gff
prodigal  -i $INFILE -o ${BASENAME}.gff  -f gff -a ${BASENAME}.proteins.faa -d ${BASENAME}.nucleotideo.fna -p meta


