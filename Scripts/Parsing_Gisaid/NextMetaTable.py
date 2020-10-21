#!/usr/bin/env python3

import argparse

'''                                              Description
    This script uses high coverage sequences (SeqH) and completed sequences from(SeqC) downloaded from Gisaid plus the 
    nextmeta.tsv also downloaded from Gisaid to reconstructed nextmeta file, but only using SeqH ^ SeqC  sequences.
    (^) = logic conjunction

'''

# using argparse to improve execution of program using terminal
parser = argparse.ArgumentParser()
parser.add_argument("sequenceFile", help="File with sequences in fasta format")
parser.add_argument("metadataFile", help="File with metadata. Could have identifiers not present in sequenceFileincrease")
parser.add_argument("outputMetadataFile",  help="Name of new file where new metadataba table will be written")
args = parser.parse_args()

# Opening files safely
try:
    seq_handler = open(args.sequenceFile, 'r')
except:
    print('File', args.sequenceFile,  'cannot be open')
    exit()

try:
    metadata_handler = open(args.metadataFile, 'r')
except:
    print('File', args.metadataFile, 'cannot be open')
    exit()

try:
    rebuild_metadata_handler = open(args.outputMetadataFile, 'w')
except:
    print('File', args.outputMetadataFile, 'cannot be open in writting mode')
    exit()

# Making list to store adds
ids = list()
for line in seq_handler:
    line = line.strip()
    if not line.startswith('>hCoV'):
        continue
    seq_id = line.split('|')[1]
    if seq_id.startswith('EPI_ISL_'):
        ids.append(seq_id)

# Creating new table to store high coverage and complete sequences
CompCoverageTable = rebuild_metadata_handler

# Making a title de um jeito bem safado
for line in metadata_handler:
    title = line
    break
CompCoverageTable.write(title)

# Creating a Hash to store identifiers
next_table = {}
for line in metadata_handler:
    string = line
    line = line.split()
    next_table[line[2]] = string

# Saving highcoverage and complete sequences in the table
matches = 0
keys_out = list()
sq_handler = open("missing_sequences.txt", "w")

for epi in ids:
    if epi in next_table:
        string = next_table[epi]
        matches += 1
        CompCoverageTable.write(string)
    else:
        keys_out.append(epi)
        sq_handler.write(epi+"\n")
# Printing run parameters
print(f'Número de matches: {matches}, IDs Totais: {len(ids)}, número de chaves fora chaves_out: {len(keys_out)}')
