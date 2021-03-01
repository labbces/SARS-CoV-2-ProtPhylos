import argparse
import os
import re
import pandas as pd
import seqlogo

# Using argparse to make easier use this script using command line
parser = argparse.ArgumentParser()
parser.add_argument("BordasPath", help="path to borders files")
parser.add_argument("Matrix",
                    help="file to save matrix")
parser.add_argument('SeqStorage', help='file to save sequences')
parser.add_argument('SeqLogo', help='file to save Sequence logo plot')
parser.add_argument('SeqType', help='aa ou nt')
parser.add_argument('NumBorder', type=int, help='Number: 1 to 15')
args = parser.parse_args()

# dictionaries to check the length of sequences and to complete the matrix
lentype = {'nt': 24,
           'aa': 8
           }
alphatype = {'nt': ['A', 'C', 'G', 'T'],
             'aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
             }

# Find and open the borders files and save borders in seqdict
patternString = '.+\.bordas\.' + str(args.NumBorder) + '_\d+_\d+\.' + args.SeqType + '\.fa$'
pattern = re.compile(patternString)
bordasFiles = os.scandir(args.BordasPath)
seqdict = {}
for file in bordasFiles:
    match = re.search(pattern, file.name)
    if match:
        filepath = args.BordasPath + '/' + file.name
        with open(filepath, "r") as file_object:
            for line in file_object:
                line = line.rstrip()
                match2 = re.search(r'^>', line)
                if not match2:
                    if len(line.upper()) == lentype[args.SeqType]:
                        seqdict[line.upper()] = 1
                    else:
                        print('The sequence size is not as expected', file.name)

# save sequences in a file and add the residues in each position at the posdict dictionary
posdict = {}
with open(args.SeqStorage, 'w') as SeqStorage:
    for seq in seqdict.keys():
        SeqStorage.write(f'{seq} \n')
        for e in range(len(seq)):
            if seq[e] in alphatype[args.SeqType]:
                if e not in posdict.keys():
                    posdict[e] = {}
                    if seq[e] not in posdict[e]:
                        posdict[e][seq[e]] = 1
                else:
                    if seq[e] not in posdict[e]:
                        posdict[e][seq[e]] = 1
                    else:
                        posdict[e][seq[e]] = posdict[e][seq[e]] + 1

# pre-matrix mold
matrix = {'nt': {'A': [],
                 'C': [],
                 'G': [],
                 'T': []
                 },
          'aa': {'A': [], 'C': [], 'D': [], 'E': [], 'F': [], 'G': [], 'H': [], 'I': [], 'K': [], 'L': [],
                 'M': [], 'N': [], 'P': [], 'Q': [], 'R': [], 'S': [], 'T': [], 'V': [], 'W': [], 'Y': []
                 }
          }

# Write matrix in to a file and build the pre-matrix dictionary
with open(args.Matrix, 'w') as outputfile:
    header = '\t'.join(sorted(alphatype[args.SeqType]))
    outputfile.write(f'  \t{header:>5}\n')
    for pos in posdict.keys():
        lenpos = sum(posdict[pos].values())
        outputfile.write(f'{pos:<2}\t')
        for res in alphatype[args.SeqType]:
            if res in sorted(posdict[pos].keys()):
                freq = posdict[pos][res] / lenpos
                matrix[args.SeqType][res].append(freq)
                outputfile.write(f'{freq:.3f}\t')
            else:
                matrix[args.SeqType][res].append(0)
                outputfile.write(f'0.000\t')
        outputfile.write('\n')

# build matrix used to make sequence logo
DataFrame = pd.DataFrame(matrix[args.SeqType])

# define alphabet type based in sequence type selected
if args.SeqType == 'nt':
    alphabet_type = 'DNA'
elif args.SeqType == 'aa':
    alphabet_type = 'AA'

# Build and save sequence logo
ppm = seqlogo.Ppm(DataFrame, alphabet_type=alphabet_type)
icFalseOut = args.SeqLogo + '.icFalse.svg'
icTrueOut = args.SeqLogo + '.icTrue.svg'

seqlogo.seqlogo(ppm, ic_scale=False, format='svg', size='medium', filename=icFalseOut)
seqlogo.seqlogo(ppm, ic_scale=True, format='svg', size='medium', filename=icTrueOut)
