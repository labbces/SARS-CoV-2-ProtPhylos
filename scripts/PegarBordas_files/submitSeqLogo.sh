#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -t 1
#$ -S /bin/bash

module load Python/3.7.2
module load pdf2svg

region=World

for type in aa nt
do
    for DatasetType in Redundant NonRedundant
    do
        for border in {1..15}
        do
            python3 Script_SequenceLogoCompleto_SemData.py ../../bordasFiles/${region}_Bordas/ Matrix_${DatasetType}_${region}_${type}_${border} Borders_${DatasetType}_${region}_${type}_${border} SequenceLogo_${DatasetType}_${region}_${type}_${border} ${DatasetType} ${type} $border
        done
    done
done
