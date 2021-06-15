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
GeneE_list = []
ProteinE_list = []

# file to save invalid files
invalid_files = open(f'invalid_files_E_{region}.txt', 'w')

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
                    if feature.type == "gene" and feature.qualifiers["gene"][0] == "E":
                        position_e = feature.location
                        name = record.id.split(".")[1]
                        e_gene = position_e.extract(record.seq)
                        e_protein = position_e.extract(record.seq).translate()
                        gene_E = SeqRecord(seq=e_gene, id="gene_E_" + name,
                                           description="envelope_protein-gene")
                        protein_E = SeqRecord(seq=e_protein, id="protein_E_" + name,
                                              description="envelope_protein-protein")
                        GeneE_list.append(gene_E)
                        ProteinE_list.append(protein_E)
                        element = feature.qualifiers["gene"][0]
                        if element == "E":
                            if element not in count1:
                                count1[element] = 1
                            else:
                                count1[element] += 1

        except:
            invalid_files.write(file)

# Printing information about how many sequences actually work
print(count1, count2)

# Writing the multifasta files with the sequences
SeqIO.write(GeneE_list, f'{output_path}/{region}_gene_e.fasta', "fasta")
SeqIO.write(ProteinE_list, f'{output_path}/{region}_protein_e.fasta', "fasta")
