import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Using argparse to improve program execution from terminal
parser = argparse.ArgumentParser()
parser.add_argument("input_Path", help="Path where the files with embl annotation are")
parser.add_argument("output_Path", help="Path where SeqI0.write will write the files")
parser.add_argument("Region", help="Region where the sequences come from")
args = parser.parse_args()

# Using path and region that was given how argument
input_path = args.input_Path
embl_path = os.scandir(input_path)
output_path = args.output_Path
region = args.Region

# List to store record from the sequences
GeneS_list = []
ProteinS_list = []

# file to save invalid files
invalid_files = open(f'invalid_files_S_{region}.txt', 'w')

# creating variables to evaluate performance
count1 = dict()
count2 = 0

# Making SeqRecords and storing in the lists
for file in embl_path:
    if file.path.endswith(".embl") and file.is_file():
        try:
            count2 += 1
            for record in SeqIO.parse(file, "embl"):
                for feature in record.features:
                    if feature.type == "gene" and feature.qualifiers["gene"][0] == "S":
                        position_s = feature.location
                        name = record.id.split(".")[1]
                        s_gene = position_s.extract(record.seq)
                        s_protein = position_s.extract(record.seq).translate()
                        gene_S = SeqRecord(seq=s_gene, id="gene_S_" + name,
                                           description="surface_glycoprotein-gene")
                        protein_S = SeqRecord(seq=s_protein, id="protein_S_" + name,
                                              description="surface_glycoprotein-protein")
                        GeneS_list.append(gene_S)
                        ProteinS_list.append(protein_S)
                        element = feature.qualifiers["gene"][0]
                        if element == "S":
                            if element not in count1:
                                count1[element] = 1
                            else:
                                count1[element] += 1

        except:
            invalid_files.write(file)

# Printing information about how many sequences actually work
print(count1, count2)

# Writing the multifasta files with the sequences
SeqIO.write(GeneS_list, f'{output_path}/{region}_gene_s.fasta', "fasta")
SeqIO.write(ProteinS_list, f'{output_path}/{region}_protein_s.fasta', "fasta")
