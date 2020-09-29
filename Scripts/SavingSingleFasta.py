from pathlib import Path
from Bio import SeqIO
import os

# Creating directories to save each sequence for region
regions = ["Oceania", "Europe", "North", "South", "Asia", "Africa"]

for region in regions:
    if not os.path.exists(region):
        os.mkdir(region)

# opening table was previously created only using highcoverage and complete seqs
table = open("HighCoverageCompleteTable.tsv", 'r')
# Creating a Hash to store identifiers
next_table = {}
for line in table:
    line = line.split()
    next_table[line[2]] = line[5]

# list to store broken sequences from Gisaid
broken = list()

# effort to save each sequence in the corresponding regions
for seq_record in SeqIO.parse("X.fasta", "fasta"):  # this fasta file was created in using copy_replace file
    try:
        y = seq_record.id.split("PIPE")[1]
        if y in next_table:
             # os.chdir(Path(f'{next_table[y]})/{seq_record.id}+.fasta')) # kind of thing I try to use
             print(next_table[y])
             SeqIO.write(seq_record, seq_record.id + "fasta", "fasta")
    except:
        print(seq_record.id)
        broken.append(seq_record.id)
print(f'número de identificadores com erro {len(broken)}')
