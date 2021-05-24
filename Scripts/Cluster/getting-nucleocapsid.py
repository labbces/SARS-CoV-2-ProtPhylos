import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Using argparse to improve program execution from terminal
parser = argparse.ArgumentParser()
parser.add_argument("Path", help="Path where the files with embl annotation are")
parser.add_argument("Region", help="Region where the sequences come from")
args = parser.parse_args()

# Using path and region that was given how argument
path = args.Path
embl_path = os.scandir(path)
region = args.Region

# List to store record from the sequences
GeneE_list = []
ProteinE_list = []

# file to save invalid files
invalid_files = open(f'invalid_files_N_{region}.txt', 'w')

# Making SeqRecords and storing in the lists
for file in embl_path:
    try:
        file = file.path
        record = SeqIO.read(file, "embl")
        position_n = record.features[35].location
        name = record.id.split(".")[1]
        n_gene = position_n.extract(record.seq)
        n_protein = position_n.extract(record.seq).translate()
        gene_N = SeqRecord(seq=n_gene, id="gene_N_" + name, description="nucleocapsid_phosphoprotein-gene")
        protein_N = SeqRecord(seq=n_protein, id="protein_N_" + name, description="nucleocapsid_phosphoprotein-protein")
        GeneE_list.append(gene_N)
        ProteinE_list.append(protein_N)
    except:
        invalid_files.write(file)

# Writing the multifasta files with the sequences
SeqIO.write(GeneE_list, f'{region}_gene_n.fasta', "fasta")
SeqIO.write(ProteinE_list, f'{region}_protein_n.fasta', "fasta")



