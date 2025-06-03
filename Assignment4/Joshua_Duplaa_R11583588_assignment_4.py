#!/usr/bin/env python3
"""
=============================================================================
Title : Joshua_Duplaa_R11583588_assignment4.py
Author : Joshua Duplaa (R11583588)
Date : 10/11/2024
Version : 3.0
Usage : python3 Joshua_Duplaa_R11583588_assignment4.py -i "in.fna" -o "out.fna" -s "BlOSUM50.mtx"
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""
####################################################################################################################################
############################################Assignment 1 code below################################################################
####################################################################################################################################
import argparse
import sys
import os
import math
import copy
import itertools


def ReadAndStoreSequences(args):
    # Validate the input file path
    if not os.path.exists(args.input):
        print(f"Error: The input file path does not exist.")
        sys.exit(1)

    with open(args.input, "r") as file:
        genome = file.readlines()
    #sequence Dictionary to store the sequence numbers and their lengths
    sequenceDict = {}
    for sequenceNum in genome:
        if(">" in sequenceNum):
            sequenceLength = countSequence(sequenceNum, 0, genome)
            sequenceDict[sequenceNum] = sequenceLength

    #Now I need to make a sorted list 
    #found cool python function for sorting the dictionary
    #https://docs.python.org/3/library/functions.html#sorted
    sorted_sequence_data = sorted(sequenceDict.items(), key=lambda x: x[1], reverse=True)
    
    return sorted_sequence_data, genome


def countSequence(sequenceNum, sequenceLength, genome):
    #looping through the relavent genome
    for i in range(genome.index(sequenceNum)+1, len(genome)):
        #check if the next sequence has not been reached. If reached return back to the main function.
        if ">" in genome[i]:
            return sequenceLength
        # Add the last sequence length if we've reached the end of the file
        elif i == len(genome) - 1:
            sequenceLength += len(genome[i].strip())
            return sequenceLength
        else:
            sequenceLength += len(genome[i].strip())  # Strip removes extra newline characters


def outSortedSeqToFile(sorted_sequence_data, genome):
    #output the sorted sequences
    sortedSequenceStrings = []

    #write sorted_sequence_data to file, but only the sequenceNum in each line. omit the length.
    with open(args.output, "w") as output_file:
        for seq in sorted_sequence_data:
            #sequence header seq[0]
            output_file.write(seq[0]) 
            genomeIndex = genome.index(seq[0])
            
            #write sequence data for each header
            for i in range(genomeIndex+1, len(genome)):
                if ">" in genome[i]:
                    break
                # Write each line of the sequence
                output_file.write(genome[i].strip() + "\n")  
                sortedSequenceStrings.append(genome[i].strip())


            #special handling for the very last sequence in the file
            if i == len(genome):
                #writing the last line of the last sequence
                output_file.write(genome[i].strip() + "\n")  
                sortedSequenceStrings.append(genome[i].strip())
            
    return sortedSequenceStrings

def grabGapscore_matrix(args):
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

    return score_matrix, gapPenalty

####################################################################################################################################
############################################Assignment 2 code below################################################################
####################################################################################################################################

def needleman_alg(sequencesToAlign, gapPenalty, score_matrix):
    #build alignMatrix
    print(sequencesToAlign)
    alignMatrix = buildAlignmentMatrix(sequencesToAlign, gapPenalty)
    traceback = buildTracebackMatrix(sequencesToAlign)
    #iterate through the matrix starting from [2, 2]
    start_row, start_col = 2, 2
    for i in range(start_row, len(alignMatrix)):
        for j in range(start_col, len(alignMatrix[i])):
            cellScore, direction = calculateMaxScore(alignMatrix, i, j, gapPenalty, score_matrix)
            alignMatrix[i][j] = cellScore
            traceback[i][j] = direction
    
    alignmentScore = alignMatrix[len(alignMatrix)-1][len(alignMatrix[0])-1]
    alignedSequences = alignSequence(traceback)

            
    return alignedSequences, alignmentScore

def buildAlignmentMatrix(sequencesToAlign, gapPenalty):
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
def calculateMaxScore(alignMatrix, rowIndex, colIndex, gapPenalty, score_matrix):
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

    #topSequence is the first sequence, topSequnce is the second sequence
    alignedSequences = [str(topSequence), str(leftSequence)]

    return alignedSequences

####################################################################################################################################
############################################Assignment 4 code below################################################################
####################################################################################################################################

def buildDistanceMatrix(sortedSequenceList, sortedSequenceStrings, score_matrix, gapPenalty, scoringType):
    distanceDict = {}
    for seq in range(len(sortedSequenceList)):
        for seqTwo in range(seq+1, len(sortedSequenceList)):
            topSequence = sortedSequenceStrings[seq]
            leftSequence = sortedSequenceStrings[seqTwo]
            sequencesToAlign = [topSequence, leftSequence]
            alignedSequences, alignmentScore = needleman_alg(sequencesToAlign, gapPenalty, score_matrix)
            distanceScore = convertToDistance(alignmentScore, topSequence, leftSequence, scoringType)
            distanceDict[seq, seqTwo] = alignmentScore, distanceScore
            #print(f"Alignment Score: {alignmentScore} ; Distance Score: {distanceScore} for sequence pair {seq+1} and {seqTwo+1} \n {alignedSequences}")
    seqCount = len(sortedSequenceStrings)
    comparisonCount = len(distanceDict)
    distanceMatrix = buildStartingDistMatrix(distanceDict, seqCount)

    print("starting distance matrix \n")
    for row in distanceMatrix:
        print(row)

    return distanceDict, distanceMatrix

def buildStartingDistMatrix(distanceDict, seqCount):
    distanceMatrix = []
    for i in range(seqCount+1):
        row = [0]*(seqCount+1)
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
    
    for key in distanceDict:
        distanceMatrix[key[0]+1][key[1]+1] = distanceDict[key][1]

    return distanceMatrix

 
def convertToDistance(alignmentScore, topSequence, leftSequence, scoringType):
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
    normalAlignmentScore = (alignmentScore-minScore)/(maxScore-minScore)
    #Convert normalized score to a disance, D = -k * ln(normalScore)
    if normalAlignmentScore <= 0:
        DistanceScore = float('inf')  # Assign a very large distance
    else:
        DistanceScore = (-1000) * math.log(normalAlignmentScore)

    return DistanceScore

#function to construct a guide tree, we only care about the mergeOrder{} dictionary though, that's the "guide tree"

def FitchMargoliash(distanceMatrix, distanceDict):
    new_distanceMatrix = copy.deepcopy(distanceMatrix)
    merged = {} #the index of this dictionary will be the order of collapse
    order = 0
    

    while len(new_distanceMatrix)>3:
        #print("\nNEW FM CLUSTERING STEP\n")
        #find the smallest distance
        min_value = min(x for row in new_distanceMatrix[1:] for x in row[1:] if x > 0)
        #print(f"min value found this iteration: {min_value}")
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
                
                    
        seqOne = new_distanceMatrix[0][rowIndex]
        seqTwo = new_distanceMatrix[0][colIndex]
        ##FIX FROM HERE ON!!
        #Merge Rows and columns of minimum and find distance from the cluster of of other sequence pairs
        #print(f"Merging sequence {seqOne} and sequence {seqTwo}")
        merged[order] = [seqOne, seqTwo]
        unmerged = new_distanceMatrix[0][1:] 
        unmerged.remove(seqOne)
        unmerged.remove(seqTwo)
        #find average between merged and unmerged
        new_distanceMatrix = remColRows(new_distanceMatrix, merged[order])

        #flatten merged to add new row and column of zeros, then fill bottom left and top right corner with merged
        newMerged = merged[order]
        newMerged = flatten(newMerged)
        #add row and column of zeros to puy newMerged on the corners
        #adding Col
        for row in range(len(new_distanceMatrix)):
            new_distanceMatrix[row] += [0]
        #adding row
        num_cols = len(new_distanceMatrix[0])
        new_distanceMatrix.append([0]*num_cols)

        new_distanceMatrix[0][num_cols-1] = newMerged
        new_distanceMatrix[num_cols-1][0] = newMerged
        
        #add averages for the proper spots
        new_distanceMatrix = findAvg(distanceDict, new_distanceMatrix)

            
        order += 1
        lastOrder = order+1
       
    newMerged = merged[order-1]
    newMerged = flatten(newMerged)
    merged[order] = [unmerged[0], newMerged]

    return new_distanceMatrix, merged

def findAvg(distanceDict, new_distanceMatrix):
    totalDist = 0 
    count = 0
    seq1 = new_distanceMatrix[0][-1]
    for row in new_distanceMatrix[1:]:
        seq2 = row[0]
        if type(seq2) == list:
            for item in seq2:
                for i in seq1:
                    if i < item:
                        totalDist += distanceDict[i-1, item-1][1]
                    elif item < i:
                        totalDist += distanceDict[item-1, i-1][1]
                    count += 1
            totalDist = totalDist/count
            #endof line index
            seq1Index = new_distanceMatrix[0].index(seq1)
            seq2Index = new_distanceMatrix[0].index(seq2)
            totalDist = new_distanceMatrix[seq2Index][seq1Index]
        else:
            for i in seq1:
                if i < seq2:
                    totalDist += distanceDict[i-1, seq2-1][1]
                elif seq2 < i:
                    totalDist += distanceDict[seq2-1, i-1][1]
                count += 1
            totalDist = totalDist/count
             #endof line index
            seq1Index = new_distanceMatrix[0].index(seq1)
            seq2Index = new_distanceMatrix[0].index(seq2)
            new_distanceMatrix[seq2Index][seq1Index] = totalDist


    return new_distanceMatrix




def flatten(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))  # Recursively flatten sublists
        else:
            result.append(item)
    return result

def remColRows(new_distanceMatrix, merged):
    #remove rows
    remMergeIndexOne = new_distanceMatrix[0].index(merged[0])
    remMergeIndexTwo = new_distanceMatrix[0].index(merged[1])
    del new_distanceMatrix[remMergeIndexOne]
    del new_distanceMatrix[remMergeIndexTwo-1]
    #remove cols
    for rowIndex in range(len(new_distanceMatrix)):
        del new_distanceMatrix[rowIndex][remMergeIndexOne]
        del new_distanceMatrix[rowIndex][remMergeIndexTwo-1]

    return new_distanceMatrix

def multipleSequenceAlignment(mergedSequences, sequenceList, gapPenalty, score_matrix):
    # Copy of original sequences
    print(sequenceList)
    alignedSequences = copy.deepcopy(sequenceList)  # Proper deep copy
    modified_sequences = copy.deepcopy(alignedSequences)  # Correct initialization
    alignmentScores = {}
    prevAligned = {}
    for item in range(len(sequenceList)):
        prevAligned[item] = set()

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
                # Perform sequence alignment
                if len(modified_sequences[i-1])>=len(modified_sequences[j-1]):
                    alignment, pairScore = msaHelper(modified_sequences[i-1], modified_sequences[j-1], gapPenalty, score_matrix)
                elif len(modified_sequences[j-1])>=len(modified_sequences[i-1]):
                    alignment, pairScore = msaHelper(modified_sequences[j-1], modified_sequences[i-1], gapPenalty, score_matrix)
                alignmentScores[order].append([i, j, pairScore, alignment])

        # Find the highest scoring alignment
        highScore = float('-inf')  # Start with negative infinity
        highScoreAlignment = None

        for item in alignmentScores[order]:
            currScore = item[2]
            if currScore > highScore:
                highScore = currScore
                highScoreAlignment = item

        print(f"High Score Alignment: {highScoreAlignment}")
        
        # Update sequences with the best alignment
        if(highScoreAlignment[0]>highScoreAlignment[1]):
            newSeqIndexOne = highScoreAlignment[1]-1
            newSeqIndexTwo = highScoreAlignment[0]-1
        else:
            newSeqIndexOne = highScoreAlignment[0]-1
            newSeqIndexTwo = highScoreAlignment[1]-1
        
        modified_sequences[newSeqIndexOne ] = highScoreAlignment[3][0]  # First aligned sequence
        modified_sequences[newSeqIndexTwo ] = highScoreAlignment[3][1]  # Second aligned sequence
        prevAligned[newSeqIndexOne].add(newSeqIndexTwo)
        prevAligned[newSeqIndexTwo].add(newSeqIndexOne)
        

        for items in prevAligned[newSeqIndexOne]:
            if items != newSeqIndexTwo:
                for index in range(len(modified_sequences[newSeqIndexOne])):
                    if modified_sequences[newSeqIndexOne][index] == "-":
                        modified_sequences[items][index] == "-"
        for items in prevAligned[newSeqIndexOne]:
            if items != newSeqIndexTwo:
                for index in range(len(modified_sequences[newSeqIndexOne])):
                    if modified_sequences[newSeqIndexOne][index] == "-":
                        modified_sequences[items][index] == "-"

        # Replace hyphens with 'X'
        modified_sequences = [seq.replace('-', 'X') for seq in modified_sequences]

        print("Matrix for next iteration:")
        print(modified_sequences)

    print(f"Aligned Sequences: {modified_sequences}")
    return modified_sequences


def msaHelper(topSequence, leftSequence, gapPenalty, score_matrix):
    #align sequences, then convert score to distance
    #swap if top is larger than left
    sequencesToAlign = topSequence, leftSequence
    AlignedSequences, AlignmentScore = needleman_alg(sequencesToAlign, gapPenalty, score_matrix)
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

def grabLabels(args):
    sequencelabels = []
    with open(args.input, "r") as file:
        genome = file.readlines()

    for row in genome:
        if row.startswith(">"):
            row = row.strip()
            sequencelabels.append(row)

    return sequencelabels

def main(args):
    #Read and put sequences in order
    sortedSequenceList, genome = ReadAndStoreSequences(args)
    sortedSequenceStrings = outSortedSeqToFile(sortedSequenceList, genome)
    score_matrix, gapPenalty = grabGapscore_matrix(args)
    scoringType = str(args.score_matrix)
    distanceDict, distanceMatrix = buildDistanceMatrix(sortedSequenceList, sortedSequenceStrings, score_matrix, gapPenalty, scoringType)
    #Do FM clustering to build guide tree from DistaneMatrix and distanceDict
    newDistanceMatrix, mergedSequences = FitchMargoliash(distanceMatrix, distanceDict)

    print(mergedSequences)
    
    #Do MSA with feng doo little 
    alignedSequencesMSA = multipleSequenceAlignment(mergedSequences, sortedSequenceStrings, gapPenalty, score_matrix)
    alignedSequencesMSA = [seq.replace('X', '-') for seq in alignedSequencesMSA]
    msaScore = sumOfPairs(alignedSequencesMSA, score_matrix)
    #grab sequence labels
    sequenceLabels = grabLabels(args)
    #now writing to file for output:: 
    with open(args.output, "w") as output_file:
        for i in range(len(alignedSequencesMSA)):
            outputString = "%s; score=%d\n%s\n" % (sequenceLabels[i], msaScore, alignedSequencesMSA[i])
            output_file.write(outputString)


    

    #mergeOrder = FitchMargoliash(distanceMatrix, alignmentDict)

#########################################

print("Assignment4 :: R#11583588")
parser = argparse.ArgumentParser(description="sorting Genome Sequences in descending order")
parser.add_argument("-i", "--input", required=True, type=str, help="File path to input.fna")
parser.add_argument("-o", "--output", required=True, type=str, help="File path for output.fna")
parser.add_argument("-s", "--score_matrix", required=True, type=str, help="File path for scoring matrix")
args = parser.parse_args()

if __name__ == "__main__":
    main(args)


#Execute task_2 function, return's clustered sequences and representative sequence list

#Execute task_3 Filter out chimeric sequeces from representative sequence list

#Execute task_4 Perform MSA of representative sequence list.