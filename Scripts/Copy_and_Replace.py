import argparse

parser = argparse.ArgumentParser()
parser.add_argument("sequenceFile", help="File with sequences in fasta format from GISAID")
parser.add_argument("newsequenceFile", help="Transformed sequence with | and / replaced")
args = parser.parse_args()

# Reading the file
file_name = args.sequenceFile
try:
    with open(file_name, 'r') as file:
        filedata = file.read()
except:
    print("File cannot be open")
    exit()

# Replacing "/" e o "|" for BARRA and PIPE respectively
newfiledata = filedata.replace('/', 'BARRA')
newfiledata = newfiledata.replace('|', 'PIPE')

# Writing new file with the changes
newseqfile = args.newsequenceFile
with open(newseqfile, 'w') as file:
    file.write(newfiledata)
    file.close()


