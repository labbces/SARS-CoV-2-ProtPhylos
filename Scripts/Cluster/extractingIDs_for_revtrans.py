from Bio import SeqIO
import argparse


# Function: Script to get IDs of sequences that represent the cluster obtained from CD-HTI

# Using argparse to use script by terminal
parser = argparse.ArgumentParser()
parser.add_argument("IDsList", help="list of ids of cluster sequences obtained by CD-HIT(proteins)")
parser.add_argument("MultiFaGenes", help="multifasta file with genes sequences")
args = parser.parse_args()

try:
    handle_id = open(args.IDsList, "r")
except:
    print("File", args.IDsList, "cannot be open")
    exit()

try:
    handle_multifa = open(args.MultiFaGenes, "r")
except:
    print("File", args.MultiFaGenes, "cannot be open")
    exit()


ids = SeqIO.parse(handle_id, "fasta")
multifa = SeqIO.parse(handle_multifa, "fasta")

gisaid_ids = dict()

for identifier in ids:
    protein_id = identifier.id.split("__")[1]
    gisaid_ids[protein_id] = 0

my_records = list()

for identifier in multifa:
    gene_id = identifier.id.split("__")[1]
    if gene_id in gisaid_ids:
        my_records.append(identifier)

SeqIO.write(my_records, "clusterGenes.fasta", "fasta")
