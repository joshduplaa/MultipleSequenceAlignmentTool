#!/usr/bin/env python3
"""
=============================================================================
Title : localAlign.py
Description : This is an implementation of  SmithWaterman to do Local pairwise alignment.
This implementation has been modified to work with the flask backend of the Pairwise Alignment web app
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 3.0
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""

def SmithWaterman(sequencesToAlign, gapPenalty, score_matrix):
    alignMatrix = BuildAlignmentMatrix(sequencesToAlign)
    traceback = BuildTracebackMatrix(sequencesToAlign)
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    cellScoreList = []
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = CalculateMaxScore(alignMatrix, i, j, gapPenalty, score_matrix)
            cellScoreList.append(cellScore)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction

    alignmentScore = max(cellScoreList)
    #find the index of the last occurrence of the cell score to determine where to start in the trace back.
    for rowIndex in range(len(alignMatrix)):  #iterate over row indices
        for colIndex in range(len(alignMatrix[rowIndex])):  #iterate over column indices
            if alignMatrix[rowIndex][colIndex] == alignmentScore:  
                #print("Encountered alignmentScore at:", rowIndex, colIndex)
                traceRow = rowIndex
                traceCol = colIndex


    alignedSequences = AlignSequence(traceback, traceRow, traceCol) #traceRow and traceCol are starting positions for the traceback
    
    #cellScore is the max value in the alignedMatrix
            
    return alignedSequences, alignmentScore



def BuildAlignmentMatrix(sequencesToAlign):
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

    #Wrinting the LEFT SEQUENCE column in matrix
    colVal = 0
    for row in traceback[2:]:
        row[0] = leftSequenceCol[colVal]
        colVal +=1

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

    #Smith Waterman condition
    if(cellScore <= 0):
        cellScore = 0
        direction = 0 #0 is the equivalent of x in smithwaterman

    return cellScore, direction

#Helper function for needleman alg to find alignment alignment between sequences
def AlignSequence(traceback, startRow, startCol):
    #Finished stepping through matrix when 0 is reached, starting point is the last occurence of the alignment score
    leftSequence = ""
    topSequence = ""
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

def Task2(sequenceList, gapPenalty, score_matrix):
    alignedSequences, alignmentScore = SmithWaterman(sequenceList, gapPenalty, score_matrix)

    return alignedSequences, alignmentScore

def main(sequenceList, gapPenalty, score_matrix):
    #Execute task 2 - aligning the sequences
    alignedSequences, alignmentScore = Task2(sequenceList, gapPenalty, score_matrix)

    return alignedSequences, alignmentScore

if __name__ == "__main__":
    main()
