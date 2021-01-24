#!/bin/bash

#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash


module load Python/3.7.2

python3 Scripts/Parsing_Gisaid/complete-coverage.py



