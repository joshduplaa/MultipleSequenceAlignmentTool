#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment2.py
Description : This is an implementation of Needleman Wunch to do global pairwise alignment.
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 3.0
Usage : python3 Joshua_Duplaa_R11583588_assignment2.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
import argparse
import sys
import os

def NeedlemanWunsch(sequencesToAlign, gapPenalty, score_matrix):
    #build alignMatrix
    alignMatrix = BuildAlignmentMatrix(sequencesToAlign, gapPenalty)
    traceback = BuildTracebackMatrix(sequencesToAlign)
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = CalculateMaxScore(alignMatrix, i, j, gapPenalty, score_matrix)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction
    
    alignmentScore = alignMatrix[len(alignMatrix)-1][len(alignMatrix[0])-1]
    alignedSequences = AlignSequence(traceback)
            
    return alignedSequences, alignmentScore

def BuildAlignmentMatrix(sequencesToAlign, gapPenalty):
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

    setRowVal = gapPenalty
    #filling in the necessary values -2,-4,-6 etc. for row 2 and beyond in alignMatrix
    for j in range(2,len(alignMatrix[0])):
        alignMatrix[1][j] = setRowVal
        setRowVal += gapPenalty    #decrement -2 everytime

    #filling in the necessary values -2,-4,-6 etc. for col 2 abd beyond in alignMatrix
    setColVal = 0
    for row in alignMatrix[1:]:
        row[1] = setColVal
        setColVal += gapPenalty

    return alignMatrix

def BuildTracebackMatrix(sequencesToAlign):
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

    #FWrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in traceback[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal +=1
    
    #filling in the necessary values L (left) for row 2 in traceback
    for j in range(2,len(traceback[0])):
        traceback[1][j] = "L"

    #filling in the necessary values U(up) for col 2 in traceback
    for row in traceback[2:]:
        row[1] = "U"

    return traceback

#Helper function for needleman alg to find alignment score.           
def CalculateMaxScore(alignMatrix, rowIndex, colIndex, gapPenalty, score_matrix):
    #find topLetter and leftLetter so we can find match/mismatch value in score_matrix
    topLetter = alignMatrix[0][colIndex]
    leftLetter = alignMatrix[rowIndex][0]

    scoreIndexTop = score_matrix[0].index(topLetter)
    scoreIndexLeft = score_matrix[0].index(leftLetter)

    #find score value at the score indexes found
    matchValue = int(score_matrix[scoreIndexTop][scoreIndexLeft])

    #Calculate score at the cell in alignMatrix
    diagVal = matchValue + alignMatrix[rowIndex-1][colIndex-1]
    upVal = alignMatrix[rowIndex-1][colIndex] + gapPenalty
    leftVal = alignMatrix[rowIndex][colIndex-1] + gapPenalty

    cellScore = max(diagVal, upVal, leftVal)
    #gap penalty found in scoring file chosen. For now just focus on the nucleatide scoring file (score_matrix)
    direction = ""
    if(cellScore == diagVal):
        direction = "D"
    elif(cellScore == upVal):
        direction = "U"
    elif(cellScore == leftVal):
        direction = "L"

    return cellScore, direction

#Helper function for needleman alg to find alignment alignment between sequences
def AlignSequence(traceback):
    #Finished stepping through matrix when arrived at traceback[1][1], starting point is the bottom left traceback[len(traceback)-1][len(traceback[0])-1]
    startRow = len(traceback)-1
    startCol = len(traceback[0])-1
    i = startRow
    j = startCol

    leftSequence = ""
    topSequence = ""

    while not (i == 1 and j == 1):
        if(traceback[i][j] == "D"):
            #Write both top and left into sequence listtopSequece
            leftSequence = traceback[i][0]+leftSequence
            topSequence = traceback[0][j]+topSequence 
            i -= 1
            j -= 1
        elif(traceback[i][j] == "U"):
            #Gap up, write left
            topSequence = "-"+topSequence 
            leftSequence = traceback[i][0]+leftSequence
            i -= 1
            
        elif(traceback[i][j] == "L"):
            #Gap Left, write top
            leftSequence = "-"+leftSequence
            topSequence = traceback[0][j]+topSequence 
            j -= 1 

    #topSequence is the first sequence, topSequnce is the second sequence
    alignedSequences = [str(topSequence), str(leftSequence)]

    return alignedSequences

def Task2(sequenceList, gapPenalty, score_matrix):
    alignedSequences, alignmentScore = NeedlemanWunsch(sequenceList, gapPenalty, score_matrix)

    return alignedSequences, alignmentScore

def main(sequenceList, gapPenalty, score_matrix):
    #Execute task 2 - aligning the sequences
    alignedSequences, alignmentScore = Task2(sequenceList, gapPenalty, score_matrix)

    return alignedSequences, alignmentScore

if __name__ == "__main__":
    main()
