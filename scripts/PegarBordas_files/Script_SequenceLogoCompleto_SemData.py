# -*- coding: utf-8 -*-

# IMPORTS
import argparse
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from logomaker import Glyph


# Using ArgParse to make easier use this script using command-line
parser = argparse.ArgumentParser()
parser.add_argument("BordasPath", help="path to borders files")
parser.add_argument("Matrix", help="file to save matrix")
parser.add_argument('SeqStorage', help='file to save sequences')
parser.add_argument('SeqLogo', help='file to save Sequence logo plot')
parser.add_argument('DataSetType', help='Redundant or NonRedundant')
parser.add_argument('SeqType', help='aa or nt')
parser.add_argument('NumBorder', type=int, help='Number: 1 to 15')
args = parser.parse_args()

# Find and Open the borders files and save borders in seqDict
lentype = {'nt': 24,
           'aa': 8
           }
seqDict = {}

patternString = '.+\.bordas\.' + str(args.NumBorder) + '_\d+_\d+\.' + args.SeqType + '\.fa$'
pattern = re.compile(patternString)
bordasFiles = os.scandir(args.BordasPath)

for file in bordasFiles:
    match = re.search(pattern, file.name)
    if match:
        filePath = args.BordasPath + '/' + file.name
        with open(filePath, "r") as file_object:
            for line in file_object:
                line = line.strip().upper()
                match2 = re.search(r'^>', line)
                if not match2:
                    if len(line) == lentype[args.SeqType]:
                        if args.DataSetType.upper().strip() == 'REDUNDANT':
                            if line in seqDict.keys():
                                seqDict[line] = seqDict[line] + 1
                            else:
                                seqDict[line] = 1
                        elif args.DataSetType.upper().strip() == 'NONREDUNDANT':
                            seqDict[line.upper()] = 1
                    else:
                        print('The sequence size is not as expected', file.name)
print(f'Dicionário das sequências:\n{seqDict}\n')


# Save sequences in a file and add the residues in each position at the seqDict dictionary
posDict = {}
alphatype = {'aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
             'nt': ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC',
                    'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT',
                    'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC',
                    'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
                    'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
             }
codons = []
ntAmount = 3

with open(args.SeqStorage, 'w') as SeqStorage:
    for seq, amount in seqDict.items():
        if args.DataSetType.upper().strip() == 'REDUNDANT':
            SeqStorage.write(f'{amount}x {seq}\n')
        elif args.DataSetType.upper().strip() == 'NONREDUNDANT':
            SeqStorage.write(f'{seq}\n')

        if args.SeqType.upper() == 'NT':
            for croppos in range(0, len(seq), ntAmount):
                codons.append(seq[croppos: croppos + ntAmount])
            for e in range(int((len(seq)) / ntAmount)):
                if codons[e] in alphatype[args.SeqType]:
                    if e not in posDict.keys():
                        posDict[e] = {codons[e]: seqDict[seq]}
                    elif e in posDict.keys():
                        if codons[e] in posDict[e].keys():
                            posDict[e][codons[e]] = posDict[e][codons[e]] + seqDict[seq]
                        else:
                            posDict[e][codons[e]] = seqDict[seq]
            codons = []

        elif args.SeqType.upper() == 'AA':
            for e in range(len(seq)):
                if seq[e] in alphatype[args.SeqType]:
                    if e not in posDict.keys():
                        posDict[e] = {}
                        if seq[e] not in posDict[e]:
                            posDict[e][seq[e]] = amount
                    else:
                        if seq[e] not in posDict[e]:
                            posDict[e][seq[e]] = amount
                        else:
                            posDict[e][seq[e]] = posDict[e][seq[e]] + 1
print(f'Dicionário das Posições:\n{posDict}\n')


# Build the matrixDict
matrixDict = {'nt': {'AAA': [], 'AAC': [], 'AAG': [], 'AAT': [], 'ACA': [], 'ACC': [], 'ACG': [], 'ACT': [], 'AGA': [],
                     'AGC': [], 'AGG': [], 'AGT': [], 'ATA': [], 'ATC': [], 'ATG': [], 'ATT': [], 'CAA': [], 'CAC': [],
                     'CAG': [], 'CAT': [], 'CCA': [], 'CCC': [], 'CCG': [], 'CCT': [], 'CGA': [], 'CGC': [], 'CGG': [],
                     'CGT': [], 'CTA': [], 'CTC': [], 'CTG': [], 'CTT': [], 'GAA': [], 'GAC': [], 'GAG': [], 'GAT': [],
                     'GCA': [], 'GCC': [], 'GCG': [], 'GCT': [], 'GGA': [], 'GGC': [], 'GGG': [], 'GGT': [], 'GTA': [],
                     'GTC': [], 'GTG': [], 'GTT': [], 'TAA': [], 'TAC': [], 'TAG': [], 'TAT': [], 'TCA': [], 'TCC': [],
                     'TCG': [], 'TCT': [], 'TGA': [], 'TGC': [], 'TGG': [], 'TGT': [], 'TTA': [], 'TTC': [], 'TTG': [],
                     'TTT': []
                     },
              'aa': {'A': [], 'C': [], 'D': [], 'E': [], 'F': [], 'G': [], 'H': [], 'I': [], 'K': [], 'L': [],
                     'M': [], 'N': [], 'P': [], 'Q': [], 'R': [], 'S': [], 'T': [], 'V': [], 'W': [], 'Y': []
                     }
          }
Codon2Symbol = {'AAA': 'A', 'AAC': 'B', 'AAG': 'C', 'AAT': 'D', 'ACA': 'E', 'ACC': 'F', 'ACG': 'G', 'ACT': 'H', 'AGA': 'I',
                    'AGC': 'J', 'AGG': 'K', 'AGT': 'L', 'ATA': 'M', 'ATC': 'N', 'ATG': 'O', 'ATT': 'P', 'CAA': 'Q', 'CAC': 'R',
                    'CAG': 'S', 'CAT': 'T', 'CCA': 'U', 'CCC': 'V', 'CCG': 'W', 'CCT': 'X', 'CGA': 'Y', 'CGC': 'Z', 'CGG': '1',
                    'CGT': '2', 'CTA': '3', 'CTC': '4', 'CTG': '5', 'CTT': '6', 'GAA': '7', 'GAC': '8', 'GAG': '9', 'GAT': ';',
                    'GCA': ':', 'GCC': '}', 'GCG': '{', 'GCT': '[', 'GGA': ']', 'GGC': '(', 'GGG': ')', 'GGT': '*', 'GTA': '&',
                    'GTC': '$', 'GTG': '%', 'GTT': '@', 'TAA': '#', 'TAC': '-', 'TAG': '=', 'TAT': '+', 'TCA': '/', 'TCC': '|',
                    'TCG': '<', 'TCT': '>', 'TGA': '?', 'TGC': '!', 'TGG': ',', 'TGT': '.', 'TTA': 'ç', 'TTC': '"', 'TTG': "'",
                    'TTT': '_'
                    }

for pos in posDict.keys():
    total = sum(posDict[pos].values())
    for res in sorted(alphatype[args.SeqType]):
        if res in posDict[pos].keys():
            freq = posDict[pos][res] / total
            matrixDict[args.SeqType][res].append(freq)
        else:
            matrixDict[args.SeqType][res].append(0)
print(f'Dicionário que formará a Matrix: \n{matrixDict}\n')

if args.SeqType.upper().strip() == 'NT':
    matrixSimb = {}
    for codon in matrixDict[args.SeqType].keys():
        matrixSimb[Codon2Symbol[codon]] = matrixDict[args.SeqType][codon]
    print(f'Dicionário que formará a Matrix com símbolos: \n{matrixSimb}\n')


# Build Matrix and save
MatraixDF = pd.DataFrame(matrixDict[args.SeqType])
MatraixProb_name = 'PROB_' + args.Matrix
MatraixDF.to_csv(sep="\t", header=True, path_or_buf=MatraixProb_name, index=True)
print(f'Matrix: \n{MatraixDF}\n')


if args.SeqType.upper().strip() == 'NT':
    matrixSymbol_name = 'PROB_' + args.Matrix + '_Symbol'
    MatraixDF = pd.DataFrame(matrixSimb)
    MatraixDF.to_csv(sep="\t", header=True, path_or_buf=matrixSymbol_name, index=True)
    print(f'Matrix de Símbolos: \n{MatraixDF}\n')



## 1 DIA -> MATRIXES PROBABILIDADE PRONTAS



# Converting probability matrix to information (bits) matrix
matrixValid = logomaker.validate_matrix(MatraixDF, matrix_type='probability', allow_nan=True)
matrixBit = logomaker.transform_matrix(matrixValid, from_type='probability', to_type='information')
matrixBit_name = "BIT_" + args.Matrix
matrixBit.to_csv(sep="\t", header=True, path_or_buf=matrixBit_name, index=True)
print(f'Matrix de Bits:\n{matrixBit}\n')

if args.SeqType.upper().strip() == 'AA':
    fig, ax = plt.subplots(figsize=[7, 3])
    # set bounding box
    ax.set_xlim([0.5, (lentype[args.SeqType]+1)])
    ax.set_ylim([0, 10])
    sequenceLogo = logomaker.Logo(matrixBit)
    # style using Logo methods
    sequenceLogo.style_xticks(anchor=0, spacing=5)
    sequenceLogo.ax.set_ylabel('information (bits)')
    sequenceLogo.ax.set_xlabel('length')
    fig.savefig('fig2.png', transparent=False)
    # FALTA SALVAR A FIGURA, MAS COMO???? ------------------------------------------------------------------

elif args.SeqType.upper().strip() == 'NT':
    dictBit = matrixBit.to_dict('Dict')
    dictBitCodon = {}
    for simb in dictBit.keys():
        for k, v in Codon2Symbol.items():
            dictBitCodon[k] = dictBit[v]
    print(f'Dict dos bits dos codons:\n{dictBitCodon}\n')

    info2Glyph = {}
    floor = ceiling = 0
    central = 2 # esse aqui eu coloquei por teste, não sei como deixar automático
    for pos in dictBitCodon[codon].keys():
        for codon in dictBitCodon.keys():
            floor = ceiling
            ceiling = (floor + dictBitCodon[codon][pos])
            if dictBitCodon[codon][pos] != 0:
                if pos not in info2Glyph.keys():
                    info2Glyph[pos] = {codon: {'bit': dictBitCodon[codon][pos], 'floor': floor, 'ceiling': ceiling, 'p': central}}
                else:
                    info2Glyph[pos][codon] = {'bit': dictBitCodon[codon][pos], 'floor': floor, 'ceiling': ceiling, 'p': central}
        central += 3
        floor = ceiling = 0
    print(f'Informações que serão usadas no glyph:\n{info2Glyph}\n')

    # removing positions and build list to fill glyphs
    PreListInfo2Glyph = list(info2Glyph.values())
    ListInfo2Glyph = []
    for codons in PreListInfo2Glyph:
        for (key, value) in codons.items():
            ListInfo2Glyph.append({
                'codon': key,
                'data': value,
            })
    print(f'Glyph final list: \n{ListInfo2Glyph}\n')

    fig, ax = plt.subplots(figsize=[7, 3])
    # set bounding box
    ax.set_xlim([0.5, (lentype[args.SeqType]+1)])
    ax.set_ylim([0, 10])

    def generate_glyph(c, p, width, ceiling, ax, floor):
        return Glyph(c=c,
                     p=p,
                     width=width,
                     ceiling=ceiling,
                     ax=ax,
                     floor=floor)
    list(
        map(
            lambda item:
            generate_glyph(
                c=item['codon'],
                p=item['data']['p'],
                width=3.0,
                ceiling=item['data']['ceiling'],
                ax=ax,
                floor=item['data']['floor'])
            ,
            ListInfo2Glyph
        )
    )

    fig.savefig('fig3.png', transparent=True)

