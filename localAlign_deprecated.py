#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment3.py
Description : This is an implementation of SmithWaterman to do local pairwise alignment.
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 3.0
Usage : python3 Joshua_Duplaa_R11583588_assignment3.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
import argparse
import sys
import os

def readAndDefineSequences(args):
    # Validate the input file path
    if not os.path.exists(args.input):
        print(f"Error: The input file path does not exist.")
        sys.exit(1)
    if not os.path.exists(args.score_matrix):
        print(f"Error: The input file path does not exist.")
        sys.exit(1)

    with open(args.input, "r") as file:
        genome = file.readlines()

    #list for sequences to align
    sequencesToAlign = ""  

    for sequenceLabel in genome:
        if(">" in sequenceLabel):
            #sequencelabel is the string before the sequence.
            sequenceString = Sequence_Store(sequenceLabel,  "", genome)            
            sequencesToAlign += sequenceString
    sequencesToAlign = sequencesToAlign.splitlines()
    return sequencesToAlign

def Sequence_Store(sequenceLabel, sequenceString, genome):
    #looping through the relavent genome file
    for i in range(genome.index(sequenceLabel)+1, len(genome)):
        #check if the next sequence has not been reached. If reached return back to the main function.
        if ">" in genome[i]:
            sequenceString += "\n"
            return sequenceString
        #add the last sequence length if we've reached the end of the file
        elif i == len(genome) - 1:
            sequenceString += genome[i]
            return sequenceString
        else:
            sequenceString += genome[i].replace("\n", "")

    return sequenceString


def buildAlignmentMatrix(sequencesToAlign):
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    alignMatrix = []

    #defining the sequences to be aligned
    topSequenceRow = sequencesToAlign[0]
    leftSequenceCol = sequencesToAlign[1]

    #they are 1 character longer already due to /n
    topSequenceLength = len(topSequenceRow)+1 #size rows 
    leftSequenceLength = len(leftSequenceCol)+1 #size of columns

    #Fill the alignMatrix with zeros
    for i in range(leftSequenceLength+1):
        row = [0] * (topSequenceLength+1)  #Create a row of zeros of length of the topsequence+1
        alignMatrix.append(row)
    
    #Writing the header TOP SEQUENCE in matrix
    #len(matrix[0]) is the length of the column
    for i in range(2,len(alignMatrix[0])):
        alignMatrix[0][i] = topSequenceRow[i-2]

    #FWrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in alignMatrix[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal += 1 


    #filling in the necessary null values 0 for row 2 and beyond in alignMatrix
    for j in range(2,len(alignMatrix[0])):
        alignMatrix[1][j] = 0

    #filling in the necessary values null values 0 for col 2 abd beyond in alignMatrix
    for row in alignMatrix[1:]:
        row[1] = 0

    return alignMatrix


def buildTracebackMatrix(sequencesToAlign):
    #Now I need to build a scoring matrix
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    traceback = []

    #defining the sequences to be aligned
    topSequenceRow = sequencesToAlign[0]
    leftSequenceCol = sequencesToAlign[1]

    #they are 1 character longer already due to /n
    topSequenceLength = len(topSequenceRow)+1 #size rows 
    leftSequenceLength = len(leftSequenceCol)+1 #size of columns

    #Fill the traceback with zeros
    for i in range(leftSequenceLength+1):
        row = [0] * (topSequenceLength+1)  #Create a row of zeros of length of the topsequence+1
        traceback.append(row)
    
    #Writing the header TOP SEQUENCE in matrix
    #len(matrix[0]) is the length of the column
    for i in range(2,len(traceback[0])):
        traceback[0][i] = topSequenceRow[i-2]

    #Wrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in traceback[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal +=1

    return traceback


def smithwaterman_alg(alignMatrix, traceback, gap_penalty, score_matrix):
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    cellScoreList = []
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = calculateMaxScore(alignMatrix, i, j, gap_penalty, score_matrix)
            cellScoreList.append(cellScore)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction

    alignmentScore = max(cellScoreList)
    #find the index of the last occurrence of the cell score to determine where to start in the trace back.
    for rowIndex in range(len(alignMatrix)):  #iterate over row indices
        for colIndex in range(len(alignMatrix[rowIndex])):  #iterate over column indices
            if alignMatrix[rowIndex][colIndex] == alignmentScore:  
                print("Encountered alignmentScore at:", rowIndex, colIndex)
                traceRow = rowIndex
                traceCol = colIndex


    alignedSequences = alignSequence(traceback, traceRow, traceCol) #traceRow and traceCol are starting positions for the traceback
    
    #cellScore is the max value in the alignedMatrix
            
    return alignedSequences, alignmentScore, alignMatrix, traceback

#Helper function for smithwaterman alg to find alignment score.           
def calculateMaxScore(alignMatrix, rowIndex, colIndex, gap_penalty, score_matrix):
    #find topLetter and leftLetter so we can find match/mismatch value in score_matrix
    topLetter = alignMatrix[0][colIndex]
    leftLetter = alignMatrix[rowIndex][0]

    scoreIndexTop = score_matrix[0].index(topLetter)
    scoreIndexLeft = score_matrix[0].index(leftLetter)

    #find score value at the score indexes found
    matchValue = int(score_matrix[scoreIndexTop][scoreIndexLeft])

    #Calculate score at the cell in alignMatrix
    diagVal = matchValue + alignMatrix[rowIndex-1][colIndex-1]
    upVal = alignMatrix[rowIndex-1][colIndex] + gap_penalty
    leftVal = alignMatrix[rowIndex][colIndex-1] + gap_penalty

    cellScore = max(diagVal, upVal, leftVal)
    #gap penalty found in scoring file chosen. For now just focus on the nucleatide scoring file (score_matrix)
    direction = ""
    if(cellScore == diagVal):
        direction = "D"
    elif(cellScore == upVal):
        direction = "U"
    elif(cellScore == leftVal):
        direction = "L"

    #Smith Waterman condition
    if(cellScore <= 0):
        cellScore = 0
        direction = 0 #0 is the equivalent of x in smithwaterman

    return cellScore, direction
#Helper function for smithwaterman alg to find alignment alignment between sequences
def alignSequence(traceback, startRow, startCol):
    #Finished stepping through matrix when 0 is reached, starting point is the last occurence of the alignment score
    leftSequence = ""
    topSequence = ""
    print(startRow)
    print(startCol)
    i = startRow
    j = startCol

    while traceback[i][j] != 0:
        if(traceback[i][j] == "D"):
            #Write both top and left into sequence listtopSequece
            leftSequence = traceback[i][0]+leftSequence
            topSequence = traceback[0][j]+topSequence 
            i -= 1
            j -= 1
        elif(traceback[i][j] == "U"):
            #Gap top, write left
            topSequence = "-"+topSequence 
            leftSequence = traceback[i][0]+leftSequence
            i -= 1
            
        elif(traceback[i][j] == "L"):
            #Gap Left, write top
            leftSequence = "-"+leftSequence
            topSequence = traceback[0][j]+topSequence 
            j -= 1 
        elif(traceback[i][j] == "0"):
            i = 1
            j = 1

    #topSequence is the first sequence, topSequnce is the second sequence
    alignedSequences = [str(topSequence), str(leftSequence)]

    return alignedSequences

    
def grabLabels(args):
    sequencelabels = []
    with open(args.input, "r") as file:
        genome = file.readlines()
    
    for row in genome:
        if row.startswith(">"):
            row = row.strip()
            sequencelabels.append(row)
    return sequencelabels

def write_matrix_to_file(matrix, filename): #helper function for writing matrices to files
    with open(filename, 'w') as file:
        for row in matrix:
            file.write(' '.join(map(str, row)) + '\n')

##############main implementation
print("Assignment3 :: R#11583588")

parser = argparse.ArgumentParser(description="sorting Genome Sequences in descending order")
parser.add_argument("-i", "--input", required=True, type=str, help="File path to input.fna")
parser.add_argument("-o", "--output", required=True, type=str, help="File path for output.fna")
parser.add_argument("-s", "--score_matrix", required=True, type=str, help="File path for scoring matrix")
args = parser.parse_args()
#storing the requested score matrix
with open(args.score_matrix, "r") as file:
    score_matrix = file.readlines()
scoreIndex = 0
for row in score_matrix:
    row = row.strip("\n")
    score_matrix[scoreIndex] = row
    scoreIndex += 1
score_matrix[0] = "0" + score_matrix[0]
score_matrix = [row.split() for row in score_matrix]

#storing gap penalty early because it is useful
gapPenalty = int(score_matrix[3][len(score_matrix[3])-1])


sequencesToAlign = readAndDefineSequences(args)
alignMatrix = buildAlignmentMatrix(sequencesToAlign)
write_matrix_to_file(alignMatrix,'align_matrix.txt')
traceback = buildTracebackMatrix(sequencesToAlign)
AlignedSequences, AlignmentScore, alignedMatrix, tracedBackMatrix = smithwaterman_alg(alignMatrix, traceback, gapPenalty, score_matrix)
#Using traceback to gap the top and left sequences
print(AlignedSequences)
print(AlignmentScore)

#grab sequence labels
sequenceLabels = grabLabels(args)
#now writing to file for output:: 
with open(args.output, "w") as output_file:
    for i in range(2):
        outputString = "%s; score=%d\n%s\n" % (sequenceLabels[i], AlignmentScore, AlignedSequences[i])
        output_file.write(outputString)

write_matrix_to_file(alignedMatrix,'ALIGNED_matrix.txt')
write_matrix_to_file(tracedBackMatrix, 'TRACED_matrix.txt')

##################end of main implementation