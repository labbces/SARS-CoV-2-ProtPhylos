#!/bin/bash

#$ -q all.q
#$ -V
#$ -cwd
#$ -pe smp 2

sed -r 's/\|EPI_.*//g' allprot1021.fasta |sed -r 's/[|/]/_/g' > allprot1021.fixnames.fasta


