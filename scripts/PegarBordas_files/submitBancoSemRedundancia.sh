#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -S /bin/bash

module load Python/3.7.2

python3 BancoSemRedundancia.py Africa_Bordas testSemRed1 nt 4
