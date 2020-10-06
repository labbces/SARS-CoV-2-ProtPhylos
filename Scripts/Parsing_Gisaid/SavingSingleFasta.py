from Bio import SeqIO
import os
import argparse

'''                                              Description
    This script uses high coverage sequences (SeqH) and completed sequences(SeqC) from Gisaid plus the reconstructed 
    metadata table, generated by the script (Next Metatable.py), to save the valid sequences in their respectively
    directories by region. 
'''

# Using argparse to set the variables by terminal easily
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


# Creating directories to save each sequence for region
regions = ["Oceania", "Europe", "North", "South", "Asia", "Africa"]
for region in regions:
    if not os.path.exists(region):
        os.mkdir(region)

# Creating a Hash to store identifiers
next_table = {}
for line in metadata_handler:
    line = line.split()
    next_table[line[2]] = line[5]

# list to store broken sequences from Gisaid
broken = list()

# effort to save each sequence in the corresponding regions
for seq_record in SeqIO.parse(seq_handler, "fasta"):  # this fasta file was created in using copy_replace file
    try:
        y = seq_record.id.split("PIPE")[1]
        if y in next_table:
             SeqIO.write(seq_record, f'{next_table[y]}/{seq_record.id}.fasta', 'fasta')
    except:
        print(seq_record.id)
        broken.append(seq_record.id)
print(f'Total of broken identifies: {len(broken)}')
