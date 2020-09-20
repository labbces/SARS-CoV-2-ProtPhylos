#!/usr/bin/env python3
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("sequenceFile", help="File with sequences in fasta format")
parser.add_argument("metadataFile", help="File with metadata. Could have identifiers not present in sequenceFileincrease")
parser.add_argument("outputMetadataFile",  help="Name of new file where new metadataba table will be written")
args = parser.parse_args()

#Global vars

seqIdsDict = dict()


#Opening files
try:
    seqFileHandler = open(args.sequenceFile, 'r')
except:
    print('File', args.sequenceFile,  'cannot be open')
    exit()

try:
    metadataFileHandler = open(args.metadataFile, 'r')
except:
    print('File',args.metadataFile, 'cannot be open')
    exit()

try:
    newMetadataFileHandler = open(args.outputMetadataFile, 'w')
except:
    print('File',args.outputMetadataFile, 'cannot be open in writting mode')
    exit()


#Process the Sequence file to get Sequences IDs in a dictionary
for line in seqFileHandler:
    line = line.strip()
    if line.startswith('>hCoV'):
#>hCoV-19/Australia/NT12/2020|EPI_ISL_426900|2020
     seqId = line.split('|')[1]
     if seqId.startswith('EPI_ISL_'):
      seqIdsDict[seqId]=1
     else:
      print("I do not understand this line, check your code!")

print('There are ', len(seqIdsDict),'sequences in file',args.sequenceFile)

####
#Process metadata table. Select the line for which there are sequence (seqIdsDict) available
countMatch=0
for line in metadataFileHandler:
 fields = line.split()
 if line.startswith("strain\tvirus"):
  newMetadataFileHandler.write(line) 
 elif fields[2] in seqIdsDict:
  countMatch += 1
  newMetadataFileHandler.write(line) 

print(countMatch,' lines in metadata file' ,args.metadataFile,'match with identifiers in', args.sequenceFile)
