import argparse
import os
import re
import pandas as pd
import seqlogo

# Using argparse to make easier use this script using command line
parser = argparse.ArgumentParser()
parser.add_argument('BordasPath', help='Caminho para os arquivos com as bordas')
parser.add_argument('MatrixStorage', help='Arquivo base para guardar a matriz')
parser.add_argument('SeqStorage', help='Arquivo para guardar todas as sequências')
parser.add_argument('SeqLogo', help='Arquivo onde será gerado o SequenceLogo')
parser.add_argument('SeqType', help='aa ou nt')
parser.add_argument('NumBorder', type=int, help='Número da borda a ser analisada. Num de 1-15')
args = parser.parse_args()

# dictionaries to check the length of sequences and to complete the matrix
lentype = {'nt': 24,
           'aa': 8
           }
alphatype = {'nt': ['A', 'C', 'G', 'T'],
             'aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
             }

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

# Find and open the borders files and add borders and frequency at the freqSeq dictionary
lsBordasFiles = os.scandir(args.BordasPath)
patternString = '.+.\.bordas\.' + str(args.NumBorder) + '_\d+_\d+\.' + args.SeqType + '\.fa$'
pattern = re.compile(patternString)
freqSeq = {}
for file in lsBordasFiles:
    match = re.search(pattern, file.name)
    if match:
        filepath = args.BordasPath + '/' + file.name
        with open(filepath, 'r') as file_object:
            for line in file_object:
                line = line.rstrip()
                match2 = re.search(r'^>', line)
                line = line.upper()
                if not match2:
                    if len(line) == lentype[args.SeqType]:
                        if line in freqSeq.keys():
                            freqSeq[line] = freqSeq[line] + 1
                        else:
                            freqSeq[line] = 1
                    else:
                        print('The sequence size is not as expected', file.name)

# To save sequences and how many times it appears in to a file
with open(args.SeqStorage, 'w') as SeqStorage:
    for key, values in freqSeq.items():
        SeqStorage.write((f'{values}x {key}\n'))

#add the residues in each position at the posdict dictionary
posDict = {}
for seq in freqSeq.keys():
    for freq in range(freqSeq[seq.upper()]):
        for pos in range(len(seq)):
            if seq[pos] in alphatype[args.SeqType]:
                if pos not in posDict.keys():
                    posDict[pos] = {}
                    if seq[pos] not in posDict[pos]:
                        posDict[pos][seq[pos]] = 1
                else:
                    if seq[pos] not in posDict[pos]:
                        posDict[pos][seq[pos]] = 1
                    else:
                        posDict[pos][seq[pos]] = posDict[pos][seq[pos]] + 1

# Write matrix in to a file and build the pre-matrix dictionary
with open(args.MatrixStorage, 'w') as outputfile:
    header = '\t'.join(sorted(alphatype[args.SeqType]))
    outputfile.write(f'  \t{header:>5}\n')
    for position in posDict.keys():
        total = sum(posDict[position].values())
        outputfile.write(f'{position:<2}\t')
        for res in alphatype[args.SeqType]: ###
            if res in sorted(posDict[position].keys()):
                freq = posDict[position][res] / total
                matrix[args.SeqType][res].append(freq)
                outputfile.write(f'{freq:.3f}\t')
            else:
                matrix[args.SeqType][res].append(0)
                outputfile.write(f'0.000\t')
        outputfile.write('\n')

# build matrix used to make sequence logo
DataFrame = pd.DataFrame(matrix[args.SeqType])
print(DataFrame)

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
