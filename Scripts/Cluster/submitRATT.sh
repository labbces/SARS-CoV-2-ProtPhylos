#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash
#$ -t 1-593

module load PAGIT/1.0

FASTA=`ls -1 dados_d_sequenciamento/Africa/fasta_files/*.fasta | head -n $SGE_TASK_ID | tail -n 1`
DIRFA=`dirname ${FASTA}`
FILEN=`basename ${FASTA}`
BASEN=${FILEN/\.fasta/.RATT}

echo ${DIRFA}
echo ${BASEN}
echo ${FASTA}

if [ -f ${DIRFA}/${BASEN}/*.final.embl ]; then
 echo "${FASTA} already processed"
else
 echo "${FASTA} processing"
 mkdir -p ${DIRFA}/${BASEN}
 cp ${FASTA} ${DIRFA}/${BASEN}
 cd ${DIRFA}/${BASEN}
 sed -i 's/\//__/g' ${FILEN}|sed 's/|/___/g'|head -n 1

 /Storage/progs/PAGIT_v1/RATT/start.ratt.sh /Storage/data2/danilo.brito/SARS-CoV-2-ProtPhylos/Embl ${FILEN} RATT Strain
 cd ../../../
fi
