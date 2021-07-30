from Bio import SeqIO
from Bio import SeqRecord
import argparse

# Using argparse to use script by terminal
parser = argparse.ArgumentParser()
parser.add_argument("MultiFaGenes", help="multifasta file with genes sequences")
args = parser.parse_args()

try:
    handle_multifa = open(args.MultiFaGenes, "r")
except:
    print("File", args.MultiFaGenes, "cannot be open")
    exit()

nc_unique = dict()

ncparser = SeqIO.parse(handle_multifa,"fasta")

for seq in ncparser:
    if seq.seq not in nc_unique.values():
        nc_unique[seq.id] = seq.seq
print(len(nc_unique.keys()))

ptr_unique = dict()

for item in nc_unique.items():
    if "*" in str(item[1].translate()[:-1]): # trimming stop codons
        continue
    if item[1].translate() not in ptr_unique.values():
        ptr_unique[item[0]] = item[1].translate()
print(len(ptr_unique.keys()))


