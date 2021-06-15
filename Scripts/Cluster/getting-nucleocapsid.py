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
GeneN_list = []
ProteinN_list = []

# file to save invalid files
invalid_files = open(f'invalid_files_N_{region}.txt', 'w')

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
                    if feature.type == "gene" and feature.qualifiers["gene"][0] == "N":
                        position_n = feature.location
                        name = record.id.split(".")[1]
                        n_gene = position_n.extract(record.seq)
                        n_protein = position_n.extract(record.seq).translate()
                        gene_N = SeqRecord(seq=n_gene, id="gene_N_" + name,
                                           description="nucleocapsid_phosphoprotein-gene")
                        protein_N = SeqRecord(seq=n_protein, id="protein_N_" + name,
                                              description="nucleocapsid_phosphoprotein-protein")
                        GeneN_list.append(gene_N)
                        ProteinN_list.append(protein_N)
                        element = feature.qualifiers["gene"][0]
                        if element == "N":
                            if element not in count1:
                                count1[element] = 1
                            else:
                                count1[element] += 1

        except:
            invalid_files.write(file)

# Printing information about how many sequences actually work
print(count1, count2)

# Writing the multifasta files with the sequences
SeqIO.write(GeneN_list, f'{output_path}/{region}_gene_n.fasta', "fasta")
SeqIO.write(ProteinN_list, f'{output_path}/{region}_protein_n.fasta', "fasta")
