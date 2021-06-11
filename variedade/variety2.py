from Bio import SeqIO
from Bio import SeqRecord
import argparse
# import re

# Function: Script to get IDs of sequences that represent the cluster obtained from CD-HTI

# Using argparse to use script by terminal
parser = argparse.ArgumentParser()
# parser.add_argument("IDsList", help="list of ids of cluster sequences obtained by CD-HIT(proteins)")
parser.add_argument("not_redundant_Ptr", help="protein without redundancy output from cd-hit run")
parser.add_argument("MultiFaGenes", help="multifasta file with genes sequences")
args = parser.parse_args()

'''try:
    handle_id = open(args.IDsList, "r")
except:
    print("File", args.IDsList, "cannot be open")
    exit()'''

try:
    handle_clstr = open(args.not_redundant_Ptr, "r")
except:
    print("File", args.not_redundant_Prt, "cannot be open")

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
    if item[1].translate() not in ptr_unique.values():
        ptr_unique[item[0]] = item[1].translate()
print(len(ptr_unique.keys()))


