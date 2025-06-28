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
import argparse
import sys
import os
import math 
#modified global aligmment script from assignment 2, returns aligned sequences and score
import util.globalAlign as globalAlign


def CalcDistanceScore(AlignmentScore, topSequence, leftSequence, scoringType):
    #first normalize the Alignment score
    #find max alignment score, which is length * match score
    #find min alignment score 
    if "nucleotide" in scoringType:
        maxMatchScore = 1
        minMissScore = -1
    elif "BLOSUM50" in scoringType:
        maxMatchScore = 13   
        minMissScore = -5
    elif "BLOSUM62" in scoringType:
        maxMatchScore = 11
        minMissScore = -4
    maxLength = max(len(topSequence), len(leftSequence)) #note, it's length of unaligned sequences
    maxScore = maxLength*maxMatchScore
    minScore = maxLength*minMissScore
    normalAlignmentScore = (AlignmentScore-minScore)/(maxScore-minScore)
    #Convert normalized score to a disance (evolutionary distance), D = -k * ln(normalScore)
    if normalAlignmentScore <= 0:
        distanceScore = float('inf')  # Assign a very large distance
    else:
        distanceScore = (-1000) * math.log(normalAlignmentScore)

    return distanceScore

def findAvgDistance(closestSequences, compareTo, distanceMatrix, distanceMatrixCopy):    
    count = 0
    sum = 0
    if type(compareTo) == list:
        for thing in compareTo:
            for item in closestSequences:
                if item<thing:
                    #add distance of the sequences being compared (thing-to-item)
                    sum += distanceMatrix[item][thing]
                else:
                    sum += distanceMatrix[thing][item]
                count += 1

    else:
        for item in closestSequences:
            if item<compareTo:
                sum += distanceMatrix[item][compareTo]
            else:
                sum += distanceMatrix[compareTo][item]
            count += 1
        

    avgDistance = sum/count

    return avgDistance

def findClosestSeqs(distanceMatrixCopy):
    newDistanceMatrix = [row[:] for row in distanceMatrixCopy[1:]]
    for row in newDistanceMatrix:   
        del row[0]
    #replace 0s with none so they aren't caught in finding closest values
    newDistanceMatrix = [[999999 if x == 0 else x for x in row] for row in newDistanceMatrix]
    
    #step 1, find the closest sequences (smallest value in the matrix), set first occurence as closest sequences
    ##TODO: optomization idea, switch to lower triangular, start from bottom
    min_value = min(min(sublist) for sublist in newDistanceMatrix)
    for row in range(len(newDistanceMatrix)):
        for col in range(len(newDistanceMatrix[row])):
            if newDistanceMatrix[row][col] == min_value:
                closestSequences = [row+1, col+1]
                break
    #return closest sequences indexes
    
    return closestSequences

def Task_2(sequenceList, gapPenalty, score_matrix, scoringType):
    #build the n+1 x n+1 empty distance matrix(List of lists), build first row with column labels (1 for seq1, 2 for seq2, etc.)
    distanceMatrix = [list(range(len(sequenceList)+1))]
    
    #create rest of rows with first item as row label
    for i in range(1,len(sequenceList)+1):
        row = [0]*(len(sequenceList)+1)
        row[0] = i
        distanceMatrix.append(row)

    #All to all alignment, store sequences and score in distance matrix(List of lists)
    for i in range(len(sequenceList)):
        for sequence in sequenceList[i+1:len(sequenceList)]:
            seq1Label = i+1
            seq2label = sequenceList.index(sequence)+1
            sequencePair = [sequenceList[i], sequence] #sequenceOne = sequenceList[i], sequenceTwo = sequence
            alignedSequences, alignmentScore = globalAlign.main(sequencePair, gapPenalty, score_matrix)
            #Calculate distance score by normalizing the alignment score and converting the normalized score to a distance
            distanceScore = CalcDistanceScore(alignmentScore, sequenceList[i], sequence, scoringType)
            #update distanceMatrix with distance score between seq1Label and seq2Label
            distanceMatrix[seq1Label][seq2label] = distanceScore

    return distanceMatrix

def Task_3(distanceMatrix):
    #create copy of distanceMatrix without the labels
    distanceMatrixCopy = [row[:] for row in distanceMatrix]
    #distance matrix takes from copy of matrix removing the labels. It's the matrix I use to find the smallest distance
    collapseOrder = []
    while(len(distanceMatrixCopy[0])>2):

        closestSequences = findClosestSeqs(distanceMatrixCopy)  #closest sequence indexes

        #Index of merged sequences
        seqToDelete1 = distanceMatrixCopy[0][closestSequences[0]]
        seqToDelete2 = distanceMatrixCopy[0][closestSequences[1]]

        indexToDelete1 = distanceMatrixCopy[0].index(seqToDelete1)
        indexToDelete2 = distanceMatrixCopy[0].index(seqToDelete2)

        #remove cols and rows for the closest sequences, then add back the merged column and row
        for row in distanceMatrixCopy:
            del row[indexToDelete1]
            if(indexToDelete1<indexToDelete2): #if statement added for edge case
                del row[indexToDelete2-1]
            else:
                del row[indexToDelete2]
        del distanceMatrixCopy[indexToDelete1]
        if(indexToDelete1<indexToDelete2):
            del distanceMatrixCopy[indexToDelete2-1]
        else:
            del distanceMatrixCopy[indexToDelete2]

        #add back merged column and row with average distances 
        ##TODO: add condition for merging clusters together, append that.
        if isinstance(seqToDelete1, list) and isinstance(seqToDelete2, list):
            merged = seqToDelete1 + seqToDelete2
        elif isinstance(seqToDelete1, list) and isinstance(seqToDelete2, int):
            merged = seqToDelete1 + [seqToDelete2]
        elif isinstance(seqToDelete1, int) and isinstance(seqToDelete2, list):
            merged = [seqToDelete1] + seqToDelete2
        elif isinstance(seqToDelete1, int) and isinstance(seqToDelete2, int):
            merged = [seqToDelete1, seqToDelete2]

         
        distanceMatrixCopy[0].append(merged) #added col with merged header
        count = 0
        
        for row in distanceMatrixCopy:
            row[0] = distanceMatrixCopy[0][count]
            count += 1
        #adding average distances into distancematrix copy
        for row in distanceMatrixCopy[1:]:
            avgDistance = findAvgDistance(merged, row[0], distanceMatrix, distanceMatrixCopy)
            row.append(avgDistance)
        #merged row
        ##TODO: need to handle merging two lists (clusters) merging
        lastRow = [merged]
        collapseOrder.append([seqToDelete1, seqToDelete2])
        for i in range(len(distanceMatrixCopy[0])-1):
            lastRow.append(0)
        distanceMatrixCopy.append(lastRow)

    return collapseOrder

def replaceGapsWithX(sequence):
    sequence = sequence.replace("-", "X")
    return sequence


def seqToList(seq, merged2, msa, gapPenalty, score_matrix):
    seenAlignments = []
    sequenceIndex1 = seq-1
    #Prior to alignment, convert all gaps to an ‘X’.
    msa[sequenceIndex1] = replaceGapsWithX(msa[sequenceIndex1])
    for seq2 in merged2:
        sequenceIndex2 = seq2-1
        #Prior to alignment, convert all gaps to an ‘X’.
        msa[sequenceIndex2] = replaceGapsWithX(msa[sequenceIndex2])
        sequencePair = [msa[sequenceIndex1], msa[sequenceIndex2]]
        alignedSequences, alignmentScore = globalAlign.main(sequencePair, gapPenalty, score_matrix)
        seenAlignments.append([alignedSequences, alignmentScore, sequenceIndex1, sequenceIndex2])       #append alignment to seenAlignments
    #Choose the pair w/ the max alignment score to replace in sequence list
    scores = []
    for item in seenAlignments:
        scores.append(item[1])
    maxScore = max(scores)
    bestAlignmentIndex = scores.index(maxScore)
    bestAlignment = seenAlignments[bestAlignmentIndex][0]

    #add new gaps from msa[seenAlignments[bestAlignmentIndex][3]] to the rest of merged 2 
    for letterIndex in range(len(bestAlignment[1])):
        if bestAlignment[1][letterIndex] == "-":
            for item in merged2:
                if(item-1 != seenAlignments[bestAlignmentIndex][3]):
                    seqIndex = item-1
                    msa[seqIndex] = msa[seqIndex][:letterIndex] + "-" + msa[seqIndex][letterIndex:]

    msa[sequenceIndex1], msa[seenAlignments[bestAlignmentIndex][3]] = bestAlignment[0], bestAlignment[1] #Updates MSA with aligned sequences

    return msa


def Task_4(collapseOrder, sequenceList, gapPenalty, score_matrix):
    #final MSA starts as a copy of sequenceList
    msa = sequenceList[:]

    for item in collapseOrder:
        merged1, merged2 = item[0], item[1]

        #Sequence to Sequence Merge
        if type(merged1) == int and type(merged2) == int:
            sequenceIndex1, sequenceIndex2 = merged1-1, merged2-1
            #Prior to alignment, convert all gaps to an ‘X’.
            msa[sequenceIndex1] = replaceGapsWithX(msa[sequenceIndex1])
            msa[sequenceIndex2] = replaceGapsWithX(msa[sequenceIndex2])
            sequencePair = [msa[sequenceIndex1], msa[sequenceIndex2]]
            alignedSequences, alignmentScore = globalAlign.main(sequencePair, gapPenalty, score_matrix) 
            msa[sequenceIndex1], msa[sequenceIndex2] = alignedSequences[0], alignedSequences[1]         #Updates MSA with aligned sequences
        

        #Sequence to MSA 
        elif type(merged1) == int and type(merged2) == list:
            msa = seqToList(merged1, merged2, msa, gapPenalty, score_matrix)
        elif type(merged1) == list and type(merged2) == int:
            msa = seqToList(merged2, merged1, msa, gapPenalty, score_matrix)
            
        #MSA vs MSA alignemnt
        elif type(merged1) == list and type(merged2) == list:
            for seq in merged1:
                seenAlignments = []
                sequenceIndex1 = seq-1
                #Prior to alignment, convert all gaps to an ‘X’.
                msa[sequenceIndex1] = replaceGapsWithX(msa[sequenceIndex1])
                for seq2 in merged2:
                    sequenceIndex2 = seq2-1
                    #Prior to alignment, convert all gaps to an ‘X’.
                    msa[sequenceIndex2] = replaceGapsWithX(msa[sequenceIndex2])
                    sequencePair = [msa[sequenceIndex1], msa[sequenceIndex2]]
                    alignedSequences, alignmentScore = globalAlign.main(sequencePair, gapPenalty, score_matrix)
                    seenAlignments.append([alignedSequences, alignmentScore, sequenceIndex1, sequenceIndex2])       #append alignment to seenAlignments
                #Choose the pair w/ the max alignment score to replace in sequence list
                scores = []
                for item in seenAlignments:
                    scores.append(item[1])
                maxScore = max(scores)
                bestAlignmentIndex = scores.index(maxScore)
                bestAlignment = seenAlignments[bestAlignmentIndex][0]
                #add new gaps from msa[sequenceIndex1] to the same positions in the rest of merged 1
                for letterIndex in range(len(bestAlignment[0])):
                    if bestAlignment[0][letterIndex] == "-":
                        for item in merged1:
                            if(item-1 != seenAlignments[bestAlignmentIndex][2]):
                                seqIndex = item-1
                                msa[seqIndex] = msa[seqIndex][:letterIndex] + "-" + msa[seqIndex][letterIndex:]

                #add new gaps from msa[seenAlignments[bestAlignmentIndex][3]] to the rest of merged 2 
                for letterIndex in range(len(bestAlignment[1])):
                    if bestAlignment[1][letterIndex] == "-":
                        for item in merged2:
                            if(item-1 != seenAlignments[bestAlignmentIndex][3]):
                                seqIndex = item-1
                                msa[seqIndex] = msa[seqIndex][:letterIndex] + "-" + msa[seqIndex][letterIndex:]
 
                
                msa[sequenceIndex1], msa[seenAlignments[bestAlignmentIndex][3]] = bestAlignment[0], bestAlignment[1] #Updates MSA with aligned sequences

    for i in range(len(msa)):
        msa[i] = [letter.replace("X","-") for letter in msa[i]]
        msa[i] = ''.join(msa[i])

    #replace "X"s back to gaps "-"
    

    return msa    #aligned sequences in the same order as the sequence list

def Task_5(alignedSequences, score_matrix):
    columnScoreList = []
    for i in range(len(alignedSequences[0])): #column number
        ## TODO: replace match and mismatchScore with what's found in score matrix
        columnScore = 0
        for j in range(len(alignedSequences)): #loops through the ith charcter of every sequeence or "row"
            current = alignedSequences[j][i]
            if current != "-":
                for k in range(j+1, len(alignedSequences)):
                    compareTo = alignedSequences[k][i]
                    indexOne, indexTwo = score_matrix[0].index(current), score_matrix[0].index(compareTo)
                    matchScore = int(score_matrix[indexOne][indexTwo]) #could be match, mismatch, or gap
                    columnScore += matchScore
        
        columnScoreList.append(int(columnScore))
                    
    msaScore = sum(columnScoreList)

    return msaScore


##############MSA util file for Final project
def main(sequenceList, gapPenalty, score_matrix):
    """
    if args.score_matrix == "input/nucleotide.mtx":
        scoringType = "nucleotide"
    elif args.score_matrix == "input/BLOSUM50.mtx":
        scoringType = "BLOSUM50"
    elif args.score_matrix == "input/BLOSUM62.mtx":
        scoringType = "BLOSUM62"
    """
    scoringType = "BLOSUM62"

    #task 2, generate distance matrix with all v all global alignemnt
    distanceMatrix = Task_2(sequenceList, gapPenalty, score_matrix, scoringType)

    #task 3, using Fitch Margoliash construct guide tree using FM clustering, return the collapse order of the distance matrix
    collapseOrder = Task_3(distanceMatrix)

    #task 4, perform MSA given collapse Order Using feng doolittle
    alignedSequences = Task_4(collapseOrder, sequenceList, gapPenalty, score_matrix)

    #task 5, calculate score of alignedSequences (Sum of Pairs!)
    msaScore = Task_5(alignedSequences, score_matrix)


    return alignedSequences, msaScore


if __name__ == "__main__":
    main()