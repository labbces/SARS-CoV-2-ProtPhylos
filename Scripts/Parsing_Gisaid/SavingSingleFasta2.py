from Bio import SeqIO
import os
import argparse

'''                                              Description
This script uses sequences from Gisaid that were filtered  by the script = complete-coverage.py, plus the reconstructed 
metadata table, generated by the script = Next Metatable2.py, to save the valid sequences in their respectively
directories by region. 
'''


# Using argparse to set the variables in the terminal
parser = argparse.ArgumentParser()
parser.add_argument("sequencesFile", help="File with sequences in fasta format")
parser.add_argument("FineMetadataFile", help="File with metadata. Contains only high coverage and complete sequencing data")
args = parser.parse_args()


# Opening files safely
try:
    seq_handler = open(args.sequencesFile, 'r+')
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
    try:
        line = line.split("\t")
        region = line[4].split('/')[0].rstrip(" ")
        identifier = line[0].split('/')[2]
        identifier = identifier.rstrip(" ")
        identifier = identifier.replace(" ", "_")
        next_table[identifier] = region
    except:
        print("oK")


# Avoiding redundancy
files_written = dict()

for root, dirs, files in os.walk("."):
    for filename in files:
        try:
            if filename.endswith(".fasta") and "hCoV-19" in filename:
                seq_names = filename.replace("____", "/").split("/")[2]
                files_written[seq_names] = 0
        except:
            print("ERROR in avoiding redundancy", filename)

# Storing broken sequences from Gisaid
broken = list()
broken_sequences = open("broken-sequences.txt", 'w')

# Counting amount of sequences that were actually written in the run
written_sequences = 0


# Using Biopython to write the sequences in their respective directory
for seq_record in SeqIO.parse(seq_handler, "fasta"):
    try:
        label = seq_record.id.split('/')[2]
        if label in files_written.keys():
            continue
        try:
            region = next_table[label]
            region = region.replace(' ', '_')
            if label in next_table:
                if not os.path.isdir(region):
                    os.mkdir(region)
                filename = seq_record.id  # configuring the path
                filename = filename.replace("/", "____")
                filename = filename.replace("|", "__")
                pathname = f'{region}/{filename}.fasta'
                SeqIO.write(seq_record, pathname, 'fasta')
                written_sequences += 1
        except:
            #print("broken identifier", seq_record.id)
            broken_sequences.write(f'broken identifier -- {seq_record.id} \n')
    except:
        # counting missing sequences and storing their respective IDs
        broken.append(seq_record.id)
        broken_sequences.write(f'broken seq.record -- {seq_record.id}  \n')

# Process's Metrics
print(f'Total of broken identifies: {len(broken)}')
print(f'Total of written sequences in the process {written_sequences}')

# Avoiding creating empty file
broken_sequences = open("broken-sequences.txt", 'r')
if len(broken_sequences.read()) == 0:
    os.remove("broken-sequences.txt")
