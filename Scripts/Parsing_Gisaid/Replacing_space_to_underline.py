#!/usr/bin/env python3
import argparse

'''                                 Description
        Script to replace space to underline to solve problem with biopython.
        Reading line by line and writing the changes.

'''
# argparse to handle variables in terminal
parser = argparse.ArgumentParser()
parser.add_argument("sequenceFile", help="File with sequences in fasta format from GISAID")
parser.add_argument("newsequenceFile", help="Transformed sequence with | and / replaced")
args = parser.parse_args()

# Reading files safely
try:
    original_seq_handler = open(args.sequenceFile, 'r')
except:
    print("File", args.sequenceFile, "cannot be open")
    exit()

try:
    transformed_sequence = open(args.newsequenceFile, 'w')
except:
    print("File", args.newsequenceFile, "cannot be open")

# Writing new file with the changes
for line in original_seq_handler:
    if line.startswith('>'):
        line = line.replace(' ', '_')  # Replacing " " to "_" to avoid problems with Biopython processing
    else:
        transformed_sequence.write(line)
transformed_sequence.close()
