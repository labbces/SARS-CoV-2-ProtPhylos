import argparse
import os
import re

parser = argparse.ArgumentParser()
parser.add_argument("BordasPath", help="arquivos com as bordas")
parser.add_argument("SequenceLogoInputFile",
                    help="arquivo a ser usado como input para o desenvolvimento do sequence logo")
parser.add_argument('SeqType', help='aa ou nt')
parser.add_argument('NumBorder', type=int, help='Número de 1-15')
args = parser.parse_args()

bordasFiles = os.scandir(args.BordasPath)
seqdict = {}
lentype = {'nt': 24,
           'aa': 8
           }
alphatype ={'nt': ['A', 'C', 'G', 'T'],
            'aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
            }
patternString = '.+\.bordas\.'+str(args.NumBorder)+'_\d+_\d+\.'+args.SeqType+'\.fa$'
pattern = re.compile(patternString)
for file in bordasFiles:
    match = re.search(pattern, file.name)
    if match:
        filepath = args.BordasPath+'/'+file.name
        with open(filepath, "r") as file_object:
            for line in file_object:
                line = line.rstrip()
                match2 = re.search(r'^>', line)
                if not match2:
                    if len(line.upper()) == lentype[args.SeqType]:
                        seqdict[line.upper()] = 1
                    else:
                        print('O tamanho da sequência não é o esperado', file.name)

        print(file.name)
posdict ={}
for seq in seqdict.keys():
    print(seq)
    for e in range(len(seq)):
        if e not in posdict.keys():
            posdict[e] = {}
            if seq[e] not in posdict[e]:
                posdict[e][seq[e]] = 1
        else:
            if seq[e] not in posdict[e]:
                posdict[e][seq[e]] = 1
            else:
                posdict[e][seq[e]] = posdict[e][seq[e]]+1

with open(args.SequenceLogoInputFile, 'w') as outputfile:
    outputfile.write(f'  \t    A\t    C\t    G\t    T\n')
    for pos in posdict.keys():
        outputfile.write(f'{pos:<2}\t')
        for res in alphatype[args.SeqType]:
            if res in posdict[pos]:
                freq = posdict[pos][res] / len(seqdict)
                outputfile.write(f'{freq:.3f}\t')
            else:
                outputfile.write(f'0.000\t')
        outputfile.write('\n')

