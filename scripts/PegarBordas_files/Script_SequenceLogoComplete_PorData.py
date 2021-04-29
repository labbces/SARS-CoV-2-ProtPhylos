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

# Metadata Table (gisaid EPI: Collection date)
metadata = {}
with open('metadata_2021-02-09_18-08.tsv', 'r') as mtdt:
    for line in mtdt:
        line = line.strip()
        match = re.search(r'^strain', line)
        if not match:
            line = line.split()
            metadata[line[2]] = line[4]
#print(metadata)

###### Prepaing data for sequence logo

# Find and Open the borders files and save borders in dataDict
lentype = {'nt': 24,
           'aa': 8
           }
dataDict = {}

patternString = '.+\.bordas\.' + str(args.NumBorder) + '_\d+_\d+\.' + args.SeqType + '\.fa$'
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
            for line in file_object:
                line = line.strip().upper()
                match2 = re.search(r'^>', line)
                if match2:
                    pattern3 = '>RATT.+(EPI_ISL_\d+).+'
                    match3 = re.search(pattern3, line)
                    gisaid_ID = (match3.group(1))
                    GISAIDdata = metadata[gisaid_ID]
                    match4 = re.search(pattern2, GISAIDdata)
                    if match4:
                        data = GISAIDdata
                        if data not in datas:
                            datas.append(data)
                    else:
                        data = 'INCOMPLETE'
                elif not match2:
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
print(f'Dicionário das datas das sequencias:\n{dataDict}\n')

datas.sort()
print(f'\nDatas {datas} \n')

# listaDados = ['2020-12-22', '1', '2', '3', '3', '4', '3', '2021-01-03']  # Lista das datas das sequências
print(len(datas))
allDays = pd.date_range(datas[0], datas[len(datas) - 1])
print(f'Todas as datas {allDays}')  # Lista de todas as datas entre a data da sequencia mais antiga e da mais nova

daysRange = 3 # pedir para o usuário !!!!!!!!!!!
initial = 0
totalGraphs = ceil(len(allDays) / daysRange)
# print(totalGraphs)

for a in range(0, totalGraphs):
    final = initial + daysRange - 1
    print(f'Initial = {initial}')
    print(f'Final = {final}')
    # print(f'Len de allDays: {len(allDays)}')
    if final > (len(allDays) - 1):
        print(f'------------------ FINAL -------------------')
        daterange = pd.date_range(allDays[initial], allDays[len(allDays) - 1])
    else:
        daterange = pd.date_range(allDays[initial], allDays[final])
    print(daterange)
    for data in daterange:
        print(data.strftime("%Y-%m-%d"))
    initial = final + 1





