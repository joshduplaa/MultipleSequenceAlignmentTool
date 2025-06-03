#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment4.py
Description : This is an implementation Feng Doolittle for progressive MSA using NW, Fitch Margialosh, and Feng Doolittle Algorithms.
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 2.0
Usage : python3 Joshua_Duplaa_R11583588_assignment4.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
import argparse
import sys
import os
import math
import copy

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


def buildAlignmentMatrix(topSequenceRow, leftSequenceCol, gapPenalty):
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    alignMatrix = []

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


def buildTracebackMatrix(topSequenceRow, leftSequenceCol):
    #Now I need to build a scoring matrix
    #Matrix of zeroes to be built below of dimensions: sequenceLength1+1 x sequenceLenght2+1
    traceback = []

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

    #Writing the LEFT SEQUENCE column in matrix
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


def needleman_alg(alignMatrix, traceback, gap_penalty, score_matrix):
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = calculateMaxScore(alignMatrix, i, j, gap_penalty, score_matrix)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction
            
    return alignMatrix, traceback

#Helper function for needleman alg to find alignment score.           
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

    return cellScore, direction
#Helper function for needleman alg to find alignment alignment between sequences
def alignSequence(traceback):
    #Finished stepping through matrix when arrived at traceback[1][1], starting point is the bottom left traceback[len(traceback)-1][len(traceback[0])-1]
    startRow = len(traceback)-1
    startCol = len(traceback[0])-1
    i = startRow
    j = startCol

    leftSequence = ""
    topSequence = ""

    while i != 1 and j != 1:
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
    
    if(j == 1):
        while(i>1):
            topSequence = "-"+topSequence 
            leftSequence = traceback[i][0]+leftSequence
            i -= 1

    if(i == 1):
        while(j>1):
            leftSequence = "-"+leftSequence
            topSequence = traceback[0][j]+topSequence 
            j -= 1 


    
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

#All of above is Assignment 2 code, below is assignment 4

#perform AlltoALl alignment
def allToAll(sequencesToAlign, gapPenalty, score_matrix):
    #Perform All to all alignment, Convert alignment scores to distances
    distanceDict = {}
    alignmentsDict = {}
    pair = 0
    for seq in range(len(sequencesToAlign)):
        for seqTwo in range(seq+1, len(sequencesToAlign)):
            topSequence = sequencesToAlign[seq]
            leftSequence = sequencesToAlign[seqTwo]
            #align sequences, then convert score to distance
            alignMatrix = buildAlignmentMatrix(topSequence, leftSequence, gapPenalty)
            traceback = buildTracebackMatrix(topSequence, leftSequence)
            alignedMatrix, tracedBackMatrix = needleman_alg(alignMatrix, traceback, gapPenalty, score_matrix)
            AlignmentScore = alignedMatrix[len(alignedMatrix)-1][len(alignedMatrix[0])-1]
            AlignedSequences = alignSequence(traceback)
            DistanceScore = convertToDistance(AlignmentScore, topSequence, leftSequence, scoringType)
            #Store Sequence pair and Distance Score in dictionary:
            distanceDict[seq, seqTwo] = DistanceScore
            print(distanceDict)
            print(AlignedSequences)
            print(f"Alignment Score: {AlignmentScore} ; Distance Score: {DistanceScore} for sequence pair {seq+1} and {seqTwo+1}")
            alignmentsDict[pair] = [[seq, seqTwo], AlignedSequences, AlignmentScore]
            pair += 1

    return distanceDict, alignmentsDict

def convertToDistance(AlignmentScore, topSequence, leftSequence, scoringType):
    #first normalize the Alignment score
    #find max alignment score, which is length * match score
    #find min alignment score 
    if 'nucleotide' in scoringType:
        maxMatchScore = 1
        minMissScore = -1
    elif 'BLOSUM50' in scoringType:
        maxMatchScore = 13   
        minMissScore = -5
    elif 'BLOSUM62' in scoringType:
        maxMatchScore = 11
        minMissScore = -4
    maxLength = max(len(topSequence), len(leftSequence))
    maxScore = maxLength*maxMatchScore
    minScore = maxLength*minMissScore
    normalAlignmentScore = (AlignmentScore-minScore)/(maxScore-minScore)
    #Convert normalized score to a disance, D = -k * ln(normalScore)
    if normalAlignmentScore <= 0:
        DistanceScore = float('inf')  # Assign a very large distance
    else:
        DistanceScore = (-1000) * math.log(normalAlignmentScore)

    return DistanceScore

def buildDistanceMatrix(DistanceDict, sequenceCount):
    #Size of the matrix is based off the number of sequences. if the sequenceCount is N, the matrice's dimension is N*N
    #initialize the distance matrix with zeros of dimension sequenceCountxsequenceCount with an extra space for labels
    distanceMatrix = []
    for i in range(sequenceCount+1):
        row = [0]*(sequenceCount+1)
        distanceMatrix.append(row)
    
    sequenceIndex = 1
    for cell in range(1, len(distanceMatrix)):
        distanceMatrix[0][cell] = sequenceIndex
        sequenceIndex += 1

    sequenceIndex = 0
    for row in distanceMatrix:
        row[0] = sequenceIndex
        sequenceIndex += 1

    print("distanceMatrix")
    for row in distanceMatrix:
        print(row)

    #fill the matrix with distances
    for (i ,j), distance in DistanceDict.items():
        distanceMatrix[i+1][j+1] = distance

    print("distanceMatrix filled")
    print(DistanceDict)
    for row in distanceMatrix:
        print(row)

        

    return distanceMatrix

def FitchMargoliash(distanceMatrix, distanceDict):
    new_distanceMatrix = copy.deepcopy(distanceMatrix)
    mergedSequences = {} #the index of this dictionary will be the order of collapse
    order = 0

    while len(new_distanceMatrix)>3:
        print("\nNEW FM CLUSTERING STEP\n")
        #find the smallest distance
        min_value = min(x for row in new_distanceMatrix[1:] for x in row[1:] if x > 0)
        print(f"min value found this iteration: {min_value}")
        xIndex = 0
        yIndex = 0
        for i in new_distanceMatrix[1:]:
            xIndex += 1
            for j in i[1:]:
                yIndex += 1
                if j == min_value:
                    rowIndex = xIndex
                    colIndex = yIndex
            yIndex = 0
                
                    
        print(f"{rowIndex}and{colIndex}")
        seqOne = new_distanceMatrix[0][rowIndex]
        seqTwo = new_distanceMatrix[0][colIndex]
        order += 1
        ##FIX FROM HERE ON!!
        #Merge Rows and columns of minimum and find distance from the cluster of of other sequence pairs
        print(f"Merging sequence {seqOne} and sequence {seqTwo}")
        new_distanceMatrix = mergeRowsColumns(new_distanceMatrix, distanceMatrix, seqOne, seqTwo)
        print(f"AFTER merging Sequences {seqOne} and {seqTwo}::")
        mergedSequences[order] = [seqOne, seqTwo]
        new_distanceMatrix[0][-1] = mergedSequences[order]
        new_distanceMatrix[-1][0] = mergedSequences[order]
        for row in new_distanceMatrix:
            print(row)

        if type(seqTwo) == list:
            newcolhead = []
            newcolhead.append(seqOne) 
            for sequence in seqTwo:
                newcolhead.append(sequence)
        else:
            newcolhead = [seqOne, seqTwo]
        

        new_distanceMatrix[0][-1] = newcolhead
        new_distanceMatrix[-1][0] = newcolhead
        
        
    print("AFTER FM ALL CLUSTERING:")
    for row in new_distanceMatrix:
        print(row)
    
    #Add the 
    order += 1
    mergedSequences[order] = [new_distanceMatrix[0][1], new_distanceMatrix[0][2]] 

    return new_distanceMatrix, mergedSequences

def mergeRowsColumns(new_distanceMatrix, distanceMatrix, i, j):
    #make copy of distance_matrix to be merged 
    mergedMatrix = copy.deepcopy(new_distanceMatrix)
    colPop = mergedMatrix[0].index(j)
    rowPop = mergedMatrix[0].index(i)

    #remove columns i and j from each row
    for row in mergedMatrix:
        row.pop(colPop)  #remove the larger index first
        row.pop(rowPop)  #then remove the smaller index

    #remove rows i and j from the matrix
    mergedMatrix.pop(colPop)
    mergedMatrix.pop(rowPop)

    #adding a merged column with averages
    mergedColHeader = [i,j]

    mergedMatrix[0] = mergedMatrix[0].append(mergedColHeader)
    for row in mergedMatrix[1:]:
        row.append(0)
        #average distance between what hasn't been merged and what has already been merged.
    
    #replace first row with column[0] values 1:
    mergedMatrix[0] = [0]*len(mergedMatrix[1])
    mergedMatrix[0][-1] = mergedColHeader
    #put copy values from first row to first col 
    for index in range(1, len(mergedMatrix[0])-1):
        mergedMatrix[0][index] = mergedMatrix[index][0]
    #add a the col header on the bottom right    
    mergedMatrix.append([0]*len(mergedMatrix[0]))
    mergedMatrix[-1][0] = mergedColHeader
    #find list of average distances to put on the last column
    for row in mergedMatrix:
        print(row)




    mergedMatrix = grabAverageDistances(distanceMatrix, mergedMatrix, mergedColHeader)

    mergedColHeader = flatten(mergedColHeader) 
 

    return mergedMatrix



def grabAverageDistances(distanceMatrix, mergedMatrix, mergedSequences):
    nonMerged = mergedMatrix[0][1:-1]
    totalDist = 0
    count = 0
    
    for i in nonMerged:
        outer = len(mergedMatrix[0])-1
        if type(i) == list:
            for seq in i:
                totalDist, count = findTotalDistCount(totalDist, count, seq, mergedSequences, distanceMatrix)
        elif type(i)==int:
            totalDist, count  = findTotalDistCount(totalDist, count, i, mergedSequences, distanceMatrix)

        
        totalDist = totalDist/count
        mergedMatrix[mergedMatrix[0].index(i)][outer] = totalDist
        count = 0
        totalDist = 0
    """
        ##edit this to reflect above.
        
            totalDist = totalDist/count    
            mergedMatrix[mergedMatrix[0].index(i)][outer] = totalDist 
            count = 0
        """
        
    return mergedMatrix

def findTotalDistCount(totalDist, count, i, mergedSequences, distanceMatrix):
    for j in mergedSequences:
        if type(j)==list:
            for item in j:
                if item > i:
                    totalDist += distanceMatrix[item][i]
                    count += 1
                else:
                    totalDist += distanceMatrix[i][item]
                    count += 1
        else:
            if j > i:
                totalDist += distanceMatrix[i][j]
                count += 1
            else:
                totalDist += distanceMatrix[j][i]
                count += 1  

    return totalDist, count  



def multipleSequenceAlignment(mergedSequences, sequenceList, gapPenalty, score_matrix):
    # Copy of original sequences
    print(sequenceList)
    alignedSequences = copy.deepcopy(sequenceList)  # Proper deep copy
    modified_sequences = copy.deepcopy(alignedSequences)  # Correct initialization
    alignmentScores = {}

    for order in mergedSequences:
        mergeSequence_1 = mergedSequences[order][0]
        mergeSequence_2 = mergedSequences[order][1]
        print(f"Aligning {mergeSequence_1} and {mergeSequence_2}")

        # Empty list to hold alignment scores for this merge step
        alignmentScores[order] = []

        # Convert to list if they are single indices
        if isinstance(mergeSequence_1, int):
            mergeSequence_1 = [mergeSequence_1]
        if isinstance(mergeSequence_2, int):
            mergeSequence_2 = [mergeSequence_2]

        for i in mergeSequence_1:
            for j in mergeSequence_2:
                print(f"{i} & {j}")

                # Perform sequence alignment
                alignment, pairScore = msaHelper(modified_sequences[i-1], modified_sequences[j-1], gapPenalty, score_matrix)
                alignmentScores[order].append([i, j, pairScore, alignment])

        # Find the highest scoring alignment
        highScore = float('-inf')  # Start with negative infinity
        highScoreAlignment = None

        for item in alignmentScores[order]:
            currScore = item[2]
            if currScore > highScore:
                highScore = currScore
                highScoreAlignment = item

        if highScoreAlignment:
            print(f"High Score Alignment: {highScoreAlignment}")

            # Update sequences with the best alignment
            newSeqIndexOne = highScoreAlignment[0]
            newSeqIndexTwo = highScoreAlignment[1]

            modified_sequences[newSeqIndexOne - 1] = highScoreAlignment[3][0]  # First aligned sequence
            modified_sequences[newSeqIndexTwo - 1] = highScoreAlignment[3][1]  # Second aligned sequence

        # Replace hyphens with 'X'
        modified_sequences = [seq.replace('-', 'X') for seq in modified_sequences]

        print("Matrix for next iteration:")
        print(modified_sequences)

    print(f"Aligned Sequences: {alignedSequences}")
    return modified_sequences


def flatten(lst):
    return sum(([flatten(x)] if isinstance(x, list) else [x] for x in lst), [])


def msaHelper(topSequence, leftSequence):
    #align sequences, then convert score to distance
    #swap if top is larger than left

    alignMatrix = buildAlignmentMatrix(topSequence, leftSequence, gapPenalty)
    traceback = buildTracebackMatrix(topSequence, leftSequence)
    alignedMatrix, tracedBackMatrix = needleman_alg(alignMatrix, traceback, gapPenalty, score_matrix)
    AlignmentScore = alignedMatrix[len(alignedMatrix)-1][len(alignedMatrix[0])-1]
    AlignedSequences = alignSequence(traceback)
    return AlignedSequences, AlignmentScore

    

def sumOfPairs(alignedSequences, score_matrix):
    #store columns as a list of lists
    SOPcols = []
    columnSequence = ""
    for col in range(len(alignedSequences[0])):
        for row in range(len(alignedSequences)):
            columnSequence += alignedSequences[row][col]
        SOPcols.append(columnSequence)
        columnSequence = ""
    #print(SOPcols)

    scoreSum = 0
    curSum = 0
    for column in SOPcols:
        while(len(column)>1):
            for letter in column[1:]:
                topLetter = column[0]
                leftLetter = letter
                scoreIndexTop = score_matrix[0].index(topLetter)
                scoreIndexLeft = score_matrix[0].index(leftLetter)
                #finds match/mismatch values
                matchValue = int(score_matrix[scoreIndexTop][scoreIndexLeft])
                curSum += matchValue
            column = column[1:]
        #print(curSum)
        scoreSum += curSum
        curSum = 0
    #print(scoreSum)
    msa_Score = scoreSum

    return msa_Score


    ##################end of main implementation
print("Assignment4 :: R#11583588")

#User arguments taken in as global variables
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
scoringType = str(args.score_matrix)

##############main implementation

#storing gap penalty early because it is useful
gapPenalty = int(score_matrix[3][len(score_matrix[3])-1])

sequencesToAlign = readAndDefineSequences(args)     #Grab list of sequences
sequenceCount = len(sequencesToAlign)
print("SequenceCount", sequenceCount)

distanceDict, alignmentsDict = allToAll(sequencesToAlign, gapPenalty, score_matrix)

#Build distance matrix
distanceMatrix = buildDistanceMatrix(distanceDict, sequenceCount)
new_distanceMatrix, mergedSequences = FitchMargoliash(distanceMatrix, distanceDict)
print(mergedSequences)
print(alignmentsDict)
msaAlignedSequences = multipleSequenceAlignment(mergedSequences, sequencesToAlign, gapPenalty, score_matrix)
msaAlignedSequences = [seq.replace('X', '-') for seq in msaAlignedSequences]

msaScore = sumOfPairs(msaAlignedSequences, score_matrix)
#grab sequence labels
sequenceLabels = grabLabels(args)
#now writing to file for output:: 
with open(args.output, "w") as output_file:
    for i in range(len(msaAlignedSequences)):
        outputString = "%s; score=%d\n%s\n" % (sequenceLabels[i], msaScore, msaAlignedSequences[i])
        output_file.write(outputString)

