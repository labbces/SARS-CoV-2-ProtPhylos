import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Using argparse to improve program execution from terminal
parser = argparse.ArgumentParser()
parser.add_argument("input_Path", help="Path where the files with embl annotation are")
parser.add_argument("output_Path", help="Path where SeqI0.write will write the files")
parser.add_argument("Region", help="Region where the sequences come from")
parser.add_argument("protein", help="Protein to be extracted")
args = parser.parse_args()

# Using path and region that was given how argument
input_path = args.input_Path
embl_path = os.scandir(input_path)
output_path = args.output_Path
region = args.Region
protein_name = args.protein

# List to store record from the sequences
Gene_list = []
Protein_list = []

# file to save invalid files
invalid_files = open(f'invalid_files_{protein_name}_{region}.txt', 'w')

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
                    if feature.type == "gene" and feature.qualifiers["gene"][0] == f'{protein_name}':
                        position = feature.location
                        name = record.id.split(".")[1]
                        cds_gene = position.extract(record.seq)
                        if len(cds_gene) % 3 != 0:
                            invalid_files.write(record.id + "\t" + "trimmed sequence" + "\n")
                            continue
                        cds_protein = position.extract(record.seq).translate()
                        gene = SeqRecord(seq=cds_gene, id=f'gene_{protein_name}_' + name)
                        protein = SeqRecord(seq=cds_protein, id=f'protein_{protein_name}_' + name)
                        Gene_list.append(gene)
                        Protein_list.append(protein)
                        element = feature.qualifiers["gene"][0]
                        if element == f'{protein_name}':
                            if element not in count1:
                                count1[element] = 1
                            else:
                                count1[element] += 1

        except:
            if file.is_file():
                invalid_files.write(file.path)

# Printing information about how many sequences actually work
print(count1, count2)

# Writing the multifasta files with the sequences
SeqIO.write(Gene_list, f'{output_path}/{region}_gene_{protein_name}.fasta', "fasta")
SeqIO.write(Protein_list, f'{output_path}/{region}_protein_{protein_name}.fasta', "fasta")
