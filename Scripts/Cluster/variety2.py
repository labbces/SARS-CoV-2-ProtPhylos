from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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
variety_for_revtrans = list()
proteins_for_alignment = list()

ncparser = SeqIO.parse(handle_multifa, "fasta")

for sequence in ncparser:
    try:
        if "*" in sequence.translate()[:-1]:  # trimming stop codons
            continue
        if sequence.seq not in nc_unique.values():
            nc_unique[sequence.id] = sequence.seq
            gene = SeqRecord(seq=sequence.seq, id=sequence.id)
            variety_for_revtrans.append(gene)
    except:
        print("ERROR: ", sequence.id)
print(len(nc_unique.keys()))

ptr_unique = dict()

for item in nc_unique.items():
    protein = SeqRecord(seq=item[1].translate(), id=item[0])
    proteins_for_alignment.append(protein)

SeqIO.write(variety_for_revtrans, 'diversidade.fasta', 'fasta')
SeqIO.write(proteins_for_alignment, 'alinhamento.fasta', 'fasta')
