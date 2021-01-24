from Bio import SeqIO
import os
import argparse

'''                                              Description
    This script uses high coverage sequences (SeqH) and completed sequences(SeqC) from Gisaid plus the reconstructed 
    metadata table, generated by the script (Next Metatable2.py), to save the valid sequences in their respectively
    directories by region. 
'''

# Using argparse to set the variables in the terminal
parser = argparse.ArgumentParser()
parser.add_argument("sequencesFile", help="File with sequences in fasta format")
parser.add_argument("FineMetadataFile", help="File with metadata. Contains only high coverage and complete sequencing data")
args = parser.parse_args()


# Opening files safely
try:
    seq_handler = open(args.sequencesFile, 'r')
except:
    print('File', args.sequencesFile,  'cannot be open')
    exit()

try:
    metadata_handler = open(args.FineMetadataFile, 'r')
except:
    print('File', args.FineMetadataFile, 'cannot be open')
    exit()


# Creating a Hash to store identifiers
next_table = {}
for line in metadata_handler:
    line = line.split("\t")
    next_table[line[2]] = line[5]

# Storing broken sequences from Gisaid
broken = list()
broken_sequences = open("broken-sequences.txt", 'w')


for seq_record in SeqIO.parse(seq_handler, "fasta"):  # Using Biopython
    try:
        y = seq_record.id.split("|")[1]
        region = next_table[y]
        region = region.replace(' ', '_')
        if y in next_table:
             if not os.path.isdir(region):
               os.mkdir(region)
             filename = seq_record.id  # configuring the path
             filename = filename.replace("|", "____")
             filename = filename.replace("/", "____")
             pathname = f'{region}/{filename}.fasta'
             SeqIO.write(seq_record, pathname, 'fasta')
    except:
        print(seq_record.id)  # counting missing sequences and storing their respective IDs
        broken.append(seq_record.id)
        broken_sequences.write(seq_record.id + "\n")
print(f'Total of broken identifies: {len(broken)}')
