# -*- coding: utf-8 -*-

# IMPORTS
import argparse
import os
import re
import numpy as np
import pandas as pd
from math import ceil
import matplotlib.pyplot as plt
import logomaker
from logomaker import Glyph
import Bio.SeqIO as seqIO


# Using ArgParse to make easier use this script using command-line
parser = argparse.ArgumentParser()
parser.add_argument("MetadataFile", help="File with metadata. Can have identifires not used in borders files")
parser.add_argument("BordasPath", help="path to borders files")
parser.add_argument("Matrix", help="file to save matrix")
parser.add_argument('SeqStorage', help='file to save sequences')
parser.add_argument('SeqLogo', help='file to save Sequence logo plot')
parser.add_argument('DataSetType', help='Redundant or NonRedundant')
parser.add_argument('SeqType', help='aa or nt')
parser.add_argument('NumBorder', type=int, help='Number: 1 to 15')
parser.add_argument('Days', type=int, help='range in days')

args = parser.parse_args()

if args.NumBorder == 11:
    quit()

# Metadata Table (gisaid ID: Collection date)
metadata = {}
with open(args.MetadataFile, 'r') as mtdt:
    for line in mtdt:
        line = line.strip()
        match = re.search(r'^strain', line)
        if not match:
            line = line.replace(" ", "")
            line = line.replace("'", "")
            line = line.split("\t")
            metadata[line[0]] = line[15]

#print(metadata)

###### Prepaing data for sequence logo

# Find and Open the borders files and save borders in dataDict
lentype = {'nt': 24,
           'aa': 8
           }
dataDict = {}

patternString = 'Bordas_.+_' + str(args.NumBorder) + '\.' + args.SeqType + '\.fasta$'
pattern = re.compile(patternString)
bordasFiles = os.scandir(args.BordasPath)
patternString2 = '\d\d\d\d-\d\d-\d\d'
pattern2 = re.compile(patternString2)

datas = []
for file in bordasFiles:
    match = re.search(pattern, file.name)
    if match:
        filePath = args.BordasPath + '/' + file.name
        with open(filePath, "r") as file_object:
            for line in file_object: #####
                match2 = re.search(r'^>', line)
                if match2:
                    match3 = re.search('>RATT\.(.+)__(.+)__(.+).final.+', line)
                    gisaid_ID = f'hCoV-19/{match3.group(1)}/{match3.group(2)}/{match3.group(3)}'
                    GISAIDdata = metadata[gisaid_ID]
                    match4 = re.search(pattern2, GISAIDdata)
                    if match4:
                        data = GISAIDdata
                        if data not in datas:
                            datas.append(data)
                    else:
                        data = 'INCOMPLETE'
                elif not match2:
                    line = line.strip().upper()
                    if data != 'INCOMPLETE':
                        if len(line) == lentype[args.SeqType]:
                            if args.DataSetType.upper().strip() == 'REDUNDANT':
                                if data not in dataDict.keys():
                                    dataDict[data] = {line: 1}
                                else:
                                    if line in dataDict[data].keys():
                                        dataDict[data][line] = dataDict[data][line] + 1
                                    else:
                                        dataDict[data][line] = 1
                            elif args.DataSetType.upper().strip() == 'NONREDUNDANT':
                                if data not in dataDict.keys():
                                    dataDict[data] = {}
                                    dataDict[data][line] = 1
                                else:
                                    dataDict[data][line] = 1
                        else:
                            print('The sequence size is not as expected', file.name)
#print(f'Dicionário das datas das datas e das sequências:\n{dataDict}\n')
datas.sort()
#print(f'\nDatas {datas} \n')

allDays = pd.date_range(datas[0], datas[len(datas) - 1])
#print(f'Todas as datas {allDays}')  # Lista de todas as datas entre a data da sequencia mais antiga e da mais nova

daysRange = args.Days
initial = 0
totalGraphs = ceil(len(allDays) / daysRange)
#print(f'Total Graphs {totalGraphs}')

DateSeqFreq = {}
for a in range(0, totalGraphs):
    final = initial + daysRange - 1
    if final > (len(allDays) - 1):
        daterange = pd.date_range(allDays[initial], allDays[len(allDays) - 1])
    else:
        daterange = pd.date_range(allDays[initial], allDays[final])
    daterangeTitle = f'{daterange[0].strftime("%Y-%m-%d")} - {daterange[(len(daterange)-1)].strftime("%Y-%m-%d")}'

    seqDict = {}
    for fulldata in daterange:
        data = fulldata.strftime("%Y-%m-%d")
        if data in datas:
            for seq, freq in dataDict[data].items():
                if seq in seqDict.keys():
                    seqDict[seq] = seqDict[seq] + freq
                else: # if seq not in seqDict.keys():
                    seqDict[seq] = freq

    if len(seqDict) != 0:  # if seqDict is not empty
        seq2Delete = []
        for seq in seqDict.keys():
            if seqDict[seq] == 1:
                seq2Delete.append(seq)
            else:
                if args.DataSetType.upper() == 'NONREDUNDANT':
                    seqDict[seq] = 1
        for seq in seq2Delete:
            del seqDict[seq]
    DateSeqFreq[daterangeTitle] = seqDict
    initial = final + 1

for GraphName, seqDict in DateSeqFreq.items():
    if len(seqDict) != 0:  # if seqDict is not empty
        posDict = {}
        alphatype = {
            'aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
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
        #print(f'Dicionário das Posições:\n{posDict}\n')

        # Build the matrixDict
        matrixDict = {
            'nt': {'AAA': [], 'AAC': [], 'AAG': [], 'AAT': [], 'ACA': [], 'ACC': [], 'ACG': [], 'ACT': [], 'AGA': [],
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
        Codon2Symbol = {'AAA': 'A', 'AAC': 'B', 'AAG': 'C', 'AAT': 'D', 'ACA': 'E', 'ACC': 'F', 'ACG': 'G', 'ACT': 'H',
                        'AGA': 'I',
                        'AGC': 'J', 'AGG': 'K', 'AGT': 'L', 'ATA': 'M', 'ATC': 'N', 'ATG': 'O', 'ATT': 'P', 'CAA': 'Q',
                        'CAC': 'R',
                        'CAG': 'S', 'CAT': 'T', 'CCA': 'U', 'CCC': 'V', 'CCG': 'W', 'CCT': 'X', 'CGA': 'Y', 'CGC': 'Z',
                        'CGG': '1',
                        'CGT': '2', 'CTA': '3', 'CTC': '4', 'CTG': '5', 'CTT': '6', 'GAA': '7', 'GAC': '8', 'GAG': '9',
                        'GAT': ';',
                        'GCA': ':', 'GCC': '}', 'GCG': '{', 'GCT': '[', 'GGA': ']', 'GGC': '(', 'GGG': ')', 'GGT': '*',
                        'GTA': '&',
                        'GTC': '$', 'GTG': '%', 'GTT': '@', 'TAA': '#', 'TAC': '-', 'TAG': '=', 'TAT': '+', 'TCA': '/',
                        'TCC': '|',
                        'TCG': '<', 'TCT': '>', 'TGA': '?', 'TGC': '!', 'TGG': ',', 'TGT': '.', 'TTA': 'ç', 'TTC': '"',
                        'TTG': "'",
                        'TTT': '_'
                        }

        for pos in posDict.keys():
            total = sum(posDict[pos].values())
            # print(f'Total: {total}')
            for res in sorted(alphatype[args.SeqType]):
                if res in posDict[pos].keys():
                    freq = posDict[pos][res] / total
                    # print(freq)
                    matrixDict[args.SeqType][res].append(freq)
                else:
                    matrixDict[args.SeqType][res].append(0)
        #print(f'Dicionário que formará a Matrix: \n{matrixDict}\n')

        if args.SeqType.upper().strip() == 'NT':
            matrixSimb = {}
            for codon in matrixDict[args.SeqType].keys():
                matrixSimb[Codon2Symbol[codon]] = matrixDict[args.SeqType][codon]
            #print(f'Dicionário que formará a Matrix com símbolos: \n{matrixSimb}\n')

        # Build Matrix and save
        MatraixDF = pd.DataFrame(matrixDict[args.SeqType])
        MatraixProb_name = (f"{GraphName}-PROB_{args.Matrix}")
        MatraixDF.to_csv(sep="\t", header=True, path_or_buf=MatraixProb_name, index=True)
        #print(f'Matrix: \n{MatraixDF}\n')

        if args.SeqType.upper().strip() == 'NT':
            matrixSymbol_name = str(GraphName) + 'PROB_' + args.Matrix + '_Symbol'
            MatraixDF = pd.DataFrame(matrixSimb)
            MatraixDF.to_csv(sep="\t", header=True, path_or_buf=matrixSymbol_name, index=True)
            #print(f'Matrix de Símbolos: \n{MatraixDF}\n')

        # Converting probability matrix to information (bits) matrix
        matrixValid = logomaker.validate_matrix(MatraixDF, matrix_type='probability', allow_nan=True)
        matrixBit = logomaker.transform_matrix(matrixValid, from_type='probability', to_type='information')
        matrixBit_name = str(GraphName) + "BIT_" + args.Matrix
        matrixBit.to_csv(sep="\t", header=True, path_or_buf=matrixBit_name, index=True)
        # print(f'Matrix de Bits:\n{matrixBit}\n')

        # Building sequence logos
        if args.SeqType.upper().strip() == 'AA':
            fig, ax = plt.subplots(figsize=[7, 4])
            # set bounding box
            ax.set_xlim([0.5, (lentype[args.SeqType] + 1)])
            ax.set_ylim([0, 10])
            sequenceLogo = logomaker.Logo(matrixBit, ax=ax)
            # style using Logo methods
            sequenceLogo.style_xticks(anchor=0, spacing=5)
            sequenceLogo.ax.set_ylabel('information (bits)')
            sequenceLogo.ax.set_xlabel('length')
            title = f'{GraphName} - {args.SeqType} {args.DataSetType} {args.NumBorder}'
            ax.set_title(title)
            figsave_name = str(GraphName) + args.SeqLogo + '.png'
            fig.savefig(figsave_name, transparent=False)

        elif args.SeqType.upper().strip() == 'NT':
            Symbols2Codon = {}
            for symbol in matrixBit.columns:
                Codon = list(Codon2Symbol.keys())[list(Codon2Symbol.values()).index(symbol)]
                Symbols2Codon[symbol] = Codon
            #print(f'Symbols 2 codon: \n {Symbols2Codon}')
            bit_matrix_codon = matrixBit.rename(Symbols2Codon, axis='columns')
            #print(f'Bit matrix codon: \n{bit_matrix_codon}\n')
            matrixCodonBit_name = str(GraphName) + "BIT_Codon_" + args.Matrix
            bit_matrix_codon.to_csv(sep="\t", header=True, path_or_buf=matrixCodonBit_name, index=True)

            dictBit = matrixBit.to_dict('Dict')
            dictBitCodon = {}
            for simb in dictBit.keys():
                for k, v in Codon2Symbol.items():
                    dictBitCodon[k] = dictBit[v]
            #print(f'Dict dos bits dos codons:\n{dictBitCodon}\n')

            color_palett = {'GCT': 'black', 'GCC': 'black', 'GCA': 'black', 'GCG': 'black', 'TTT': 'black',
                            'TTC': 'black',
                            'ATT': 'black', 'ATC': 'black', 'ATA': 'black', 'TTA': 'black', 'TTG': 'black',
                            'CTT': 'black',
                            'CTC': 'black', 'CTA': 'black', 'CTG': 'black', 'ATG': 'black', 'CCT': 'black',
                            'CCC': 'black',
                            'CCA': 'black', 'CCG': 'black', 'TGG': 'black', 'GTT': 'black', 'GTC': 'black',
                            'GTA': 'black',
                            'GTG': 'black',
                            'TGT': 'lime', 'TGC': 'lime', 'GGT': 'lime', 'GGC': 'lime', 'GGA': 'lime', 'GGG': 'lime',
                            'TCT': 'lime', 'TCC': 'lime', 'TCA': 'lime', 'TCG': 'lime', 'AGT': 'lime', 'AGC': 'lime',
                            'ACT': 'lime', 'ACC': 'lime', 'ACA': 'lime', 'ACG': 'lime', 'TAT': 'lime', 'TAC': 'lime',
                            'GAT': 'red', 'GAC': 'red', 'GAA': 'red', 'GAG': 'red',
                            'CAT': 'blue', 'CAC': 'blue', 'AAA': 'blue', 'AAG': 'blue', 'CGT': 'blue', 'CGC': 'blue',
                            'CGA': 'blue', 'CGG': 'blue', 'AGA': 'blue', 'AGG': 'blue',
                            'AAT': 'magenta', 'AAC': 'magenta', 'CAA': 'magenta', 'CAG': 'magenta',
                            'TAA': 'peru', 'TAG': 'gold', 'TGA': 'lightcyan'}

            info2Glyph = {}
            p = 2
            for pos in range(0, int((lentype[args.SeqType]) / 3)):
                # print(f'pos: \n{pos}\n')
                codonbits = {}
                floor = ceiling = 0
                for codon in dictBitCodon.keys():
                    codonbits[codon] = dictBitCodon[codon][pos]
                codonBits_Order = sorted(codonbits.items(), key=lambda x: x[1], reverse=False)
                for TupleCodonBit in codonBits_Order:
                    trinca = TupleCodonBit[0]
                    bit = TupleCodonBit[1]
                    if bit != 0.0:
                        # print(f'TupleCodonBit: {TupleCodonBit}')
                        floor = ceiling
                        ceiling = (floor + TupleCodonBit[1])
                        # print(f'Floor: {floor}')
                        # print(f'Ceiling: {ceiling}')
                        if pos not in info2Glyph.keys():
                            info2Glyph[pos] = {trinca: {'bit': bit, 'floor': floor, 'ceiling': ceiling, 'p': p,
                                                        'color': color_palett[trinca]}}
                        else:
                            info2Glyph[pos][trinca] = {'bit': bit, 'floor': floor, 'ceiling': ceiling, 'p': p,
                                                       'color': color_palett[trinca]}
                p += 3
            #print(f'\ninfo2Glyph: \n{info2Glyph}\n')

            # removing positions and build list to fill glyphs
            PreListInfo2Glyph = list(info2Glyph.values())
            ListInfo2Glyph = []
            for codons in PreListInfo2Glyph:
                for (key, value) in codons.items():
                    ListInfo2Glyph.append({
                        'codon': key,
                        'data': value,
                    })
            #print(f'Glyph final list: \n{ListInfo2Glyph}\n')

            fig, ax = plt.subplots(figsize=[7, 4])
            # set bounding box
            ax.set_xlim([0.5, (lentype[args.SeqType] + 1)])
            ax.set_ylim([0, 8])


            def generate_glyph(c, p, width, ceiling, ax, floor, color):
                return Glyph(c=c,
                             p=p,
                             width=width,
                             ceiling=ceiling,
                             ax=ax,
                             floor=floor,
                             color=color)


            list(
                map(
                    lambda item:
                    generate_glyph(
                        c=item['codon'],
                        p=item['data']['p'],
                        width=3.0,
                        ceiling=item['data']['ceiling'],
                        ax=ax,
                        floor=item['data']['floor'],
                        color=item['data']['color'])
                    ,
                    ListInfo2Glyph
                )
            )
            title = f'{GraphName} - {args.SeqType} {args.DataSetType} {args.NumBorder}'
            ax.set_title(title)
            ax.set_ylabel('information (bits)')
            ax.set_xlabel('length')
            figsave_name = str(GraphName) + args.SeqLogo + '.png'
            fig.savefig(figsave_name, transparent=False)
    # else:
        # +(f'Date range without sequences: {daterange}')
