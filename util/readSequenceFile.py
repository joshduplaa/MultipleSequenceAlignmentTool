#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment4.py
Description : This is an implementation Feng Doolittle for progressive MSA using NW, Fitch Margialosh clustering, and sum of pairs.
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 2.0
Usage : python3 Joshua_Duplaa_R11583588_assignment4.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
import sys
import os
#modified global aligmment script from assignment 2, returns aligned sequences and score
def ReadAndStoreSequences(inputFile):
    #sequenceLabel has the title, sequenceDict has title:sequence sequenceDict[title]=sequence
    #validate the input file path
    if not os.path.exists(inputFile):
        print(f"Error: The input file path does not exist.")
        sys.exit(1)

    with open(inputFile, "r") as file:
        sequenceFile = [line.strip() for line in file]

    #dictionary for quick sequence lookup
    labelToSeq = {}
    sequenceList = []
    titleList = []
    seqIndex = -1
 
    for line in sequenceFile:
        if line[0] == ">":
            #skip line and set label To Seq
            label = line
            sequence = ""
            sequenceList.append("")
            titleList.append(line)
            seqIndex += 1
        else:
            sequence += line
            labelToSeq[label] = sequence
            sequenceList[seqIndex] = sequence
    
    seqToLabel = {value: key for key, value in labelToSeq.items()}
    #sortedSequenceList = sorted(sequenceList, key=len, reverse=True)

    return sequenceList, labelToSeq, seqToLabel, titleList

##############main implementation
def main(inputFile):
    #read in the sequences
    sequenceList, labelToSeq, seqToLabel, titleList = ReadAndStoreSequences(inputFile)
    print(sequenceList)
    return sequenceList, labelToSeq, seqToLabel, titleList

if __name__ == "__main__":
    main()