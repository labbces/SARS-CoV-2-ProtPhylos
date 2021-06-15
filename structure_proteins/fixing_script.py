from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

hash = dict()
structural = ["E", "M", "N", "S"]
count = 0

directory = r'/home/dan/sarsic/SARS-CoV-2-ProtPhylos/teste_2'
for entry in os.scandir(directory):
    if entry.path.endswith(".embl") and entry.is_file():
        count += 1
        for record in SeqIO.parse(entry, "embl"):
            for feature in record.features:
                for protein in structural:
                    if feature.type == "gene" and feature.qualifiers["gene"][0] == protein:
                        element = feature.qualifiers["gene"][0]
                        if element in ["S", "M", "N", "E"]:
                            if element not in hash:
                                hash[element] = 1
                            else:
                                hash[element] += 1

print(hash, count)

