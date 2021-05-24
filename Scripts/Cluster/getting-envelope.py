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
invalid_files = open(f'invalid_files_E_{region}.txt', 'w')

# Making SeqRecords and storing in the lists
for file in embl_path:
    try:
        file = file.path
        record = SeqIO.read(file, "embl")
        position_e = record.features[25].location
        name = record.id.split(".")[1]
        e_gene = position_e.extract(record.seq)
        e_protein = position_e.extract(record.seq).translate()
        gene_E = SeqRecord(seq=e_gene, id="gene_E_" + name, description="envelope_protein-gene")
        protein_E = SeqRecord(seq=e_protein, id="protein_E_" + name, description="envelope_protein-protein")
        GeneE_list.append(gene_E)
        ProteinE_list.append(protein_E)
    except:
        invalid_files.write(file)

# Writing the multifasta files with the sequences
SeqIO.write(GeneE_list, f'{region}_gene_e.fasta', "fasta")
SeqIO.write(ProteinE_list, f'{region}_protein_e.fasta', "fasta")



