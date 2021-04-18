from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse


# Using argparse to improve program execution from terminal
parser = argparse.ArgumentParser()
parser.add_argument("AnnotationFile", help="Annotation from RATT")
args = parser.parse_args()

# Opening file safely
try:
    handle = open(args.AnnotationFile, 'r')
except:
    print('File', args.AnnotationFile, 'cannot be open')
    exit()

# Using SeqIO in read mode
record = SeqIO.read(handle, "embl")
name = record.id.split(".")[1]

# Making SeqRecords
position_s = record.features[21].location
s_gene = position_s.extract(record.seq)
s_protein = position_s.extract(record.seq).translate()
gene_S = SeqRecord(seq=s_gene, id="gene_S_"+name, description="surface_glycoprotein-gene")
protein_S = SeqRecord(seq=s_protein, id="protein_S_"+name, description="surface_glycoprotein-protein")
position_e = record.features[25].location
e_gene = position_e.extract(record.seq)
e_protein = position_e.extract(record.seq).translate()
gene_E = SeqRecord(seq=e_gene, id="gene_E_" + name, description="envelope_protein-gene")
protein_E = SeqRecord(seq=e_protein, id="protein_E_" + name, description="envelope_protein-protein")
position_m = record.features[27].location
m_gene = position_m.extract(record.seq)
m_protein = position_m.extract(record.seq).translate()
gene_M = SeqRecord(seq=m_gene, id="gene_M_" + name, description="membrane_glycoprotein-gene")
protein_M = SeqRecord(seq=m_protein, id="protein_M_" + name, description="membrane glycoprotein-protein")
position_n = record.features[35].location
n_gene = position_n.extract(record.seq)
n_protein = position_n.extract(record.seq).translate()
gene_N = SeqRecord(seq=n_gene, id="gene_N_" + name, description="nucleocapsid_phosphoprotein-gene")
protein_N = SeqRecord(seq=n_protein, id="protein_N_" + name, description="nucleocapsid_phosphoprotein-protein")

# Writing structural gene and protein
SeqIO.write(protein_S, "protein_S_" + name, "fasta")
SeqIO.write(gene_S, "gene_S_" + name, "fasta")
SeqIO.write(protein_E, "protein_E_" + name, "fasta")
SeqIO.write(gene_E, "gene_E_" + name, "fasta")
SeqIO.write(protein_M, "protein_M_" + name, "fasta")
SeqIO.write(gene_M, "gene_M_" + name, "fasta")
SeqIO.write(protein_N, "protein_N_"+name, "fasta")
SeqIO.write(gene_N, "gene_N_"+name, "fasta")
