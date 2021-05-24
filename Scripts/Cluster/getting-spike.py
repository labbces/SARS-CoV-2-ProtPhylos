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
invalid_files = open(f'invalid_files_S_{region}.txt', 'w')

# Making SeqRecords and storing in the lists
for file in embl_path:
    try:
        file = file.path
        record = SeqIO.read(file, "embl")
        position_s = record.features[21].location
        name = record.id.split(".")[1]
        s_gene = position_s.extract(record.seq)
        s_protein = position_s.extract(record.seq).translate()
        gene_S = SeqRecord(seq=s_gene, id="gene_S_" + name, description="surface_glycoprotein-gene")
        protein_S = SeqRecord(seq=s_protein, id="protein_S_" + name, description="surface_glycoprotein-protein")
        GeneE_list.append(gene_S)
        ProteinE_list.append(protein_S)
    except:
        invalid_files.write(file)

# Writing the multifasta files with the sequences
SeqIO.write(GeneE_list, f'{region}_gene_s.fasta', "fasta")
SeqIO.write(ProteinE_list, f'{region}_protein_s.fasta', "fasta")



