import argparse

'''                                 Description
        Copy and Replace multifasta file from Gisaid. Fixing the problem to save in reason of "/" and "Pipe".
        Reading line by line and writing the changes.

'''
# argparse to handle variables in terminal
parser = argparse.ArgumentParser()
parser.add_argument("sequenceFile", help="File with sequences in fasta format from GISAID")
parser.add_argument("newsequenceFile", help="Transformed sequence with | and / replaced")
args = parser.parse_args()

# Reading files safely
try:
    original_seq_handler = open(args.sequenceFile, 'r')
except:
    print("File", args.sequenceFile, "cannot be open")
    exit()

try:
    transformed_sequence = open(args.newsequenceFile, 'w')
except:
    print("File", args.newsequenceFile, "cannot be open")

# Writing new file with the changes
for line in original_seq_handler:
    if line.startswith('>'):
        line = line.replace(' ', '_')  # Replacing " " for "_"
        transformed_sequence.write(line)
    else:
        transformed_sequence.write(line)
transformed_sequence.close()