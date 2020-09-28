import argparse

parser = argparse.ArgumentParser()
parser.add_argument("sequenceFile", help="File with sequences in fasta format")
parser.add_argument("metadataFile", help="File with metadata. Could have identifiers not present in sequenceFileincrease")
args = parser.parse_args()

# Opening files safely
try:
    seq_handler = open(args.sequenceFile, 'r')
except:
    print('Sequence file cannot be open')
    exit()

try:
    meta_handler = open(args.metadataFile, 'r')
except:
    print('Metadate file cannot be open')
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
CompCoverageTable = open('HighCoverageCompleteTable.tsv', 'w')

# Making a title
for line in meta_handler:
    title = line
    break
CompCoverageTable.write(title)

# Creating a Hash to store identifiers
next_table = {}
for line in meta_handler:
    string = line
    line = line.split()
    next_table[line[2]] = string

# Saving highcoverage and complete sequences in the table
matches = 0
chaves_out = list()
for epi in ids:
    if epi in next_table:
        string = next_table[epi]
        matches += 1
        CompCoverageTable.write(string)
    else:
        chaves_out.append(epi)

# Printing run parameters
print(f'Número de matches: {matches}, IDs Totais: {len(ids)}, número de chaves fora chaves_out: {len(chaves_out)}')