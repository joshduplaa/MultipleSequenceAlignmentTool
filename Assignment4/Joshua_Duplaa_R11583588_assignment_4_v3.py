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

import math
import copy

def buildDistanceMatrix(sorted_sequence_strings, score_matrix, gap_penalty, scoring_type):
    distance_dict = {}
    seq_count = len(sorted_sequence_strings)
    
    for i in range(seq_count):
        for j in range(i + 1, seq_count):
            top_sequence = sorted_sequence_strings[i]
            left_sequence = sorted_sequence_strings[j]
            sequences_to_align = [top_sequence, left_sequence]
            
            aligned_sequences, alignment_score = needleman_alg(sequences_to_align, gap_penalty, score_matrix)
            distance_score = convertToDistance(alignment_score, top_sequence, left_sequence, scoring_type)
            distance_dict[(i, j)] = alignment_score, distance_score
    
    distance_matrix = build_starting_dist_matrix(distance_dict, seq_count)
    return distance_dict, distance_matrix

def build_starting_dist_matrix(distance_dict, seq_count):
    distance_matrix = [[0] * (seq_count + 1) for _ in range(seq_count + 1)]
    
    for i in range(1, seq_count + 1):
        distance_matrix[0][i] = i
        distance_matrix[i][0] = i
    
    for (i, j), (_, dist) in distance_dict.items():
        distance_matrix[i + 1][j + 1] = dist
        distance_matrix[j + 1][i + 1] = dist
    
    return distance_matrix
 
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
    merged = {}  # Dictionary to store merge order
    order = 0

    while len(new_distanceMatrix) > 3:
        # Find the smallest distance
        min_value = min(x for row in new_distanceMatrix[1:] for x in row[1:] if x > 0)
        rowIndex, colIndex = None, None

        for i in range(1, len(new_distanceMatrix)):
            for j in range(1, len(new_distanceMatrix[i])):
                if new_distanceMatrix[i][j] == min_value:
                    rowIndex, colIndex = i, j
                    break
            if rowIndex is not None:
                break

        seqOne, seqTwo = new_distanceMatrix[0][rowIndex], new_distanceMatrix[0][colIndex]
        merged[order] = [seqOne, seqTwo]
        unmerged = [seq for seq in new_distanceMatrix[0][1:] if seq not in merged[order]]

        # Remove columns and rows of the merged sequences
        new_distanceMatrix = remColRows(new_distanceMatrix, merged[order])

        # Flatten merged sequences
        newMerged = flatten(merged[order])
        
        # Extend matrix with a new row and column for the merged sequence
        for row in new_distanceMatrix:
            row.append(0)
        new_distanceMatrix.append([0] * len(new_distanceMatrix[0]))
        
        new_distanceMatrix[0][-1] = newMerged
        new_distanceMatrix[-1][0] = newMerged
        
        # Compute new distances
        new_distanceMatrix = findAvg(distanceDict, new_distanceMatrix)
        
        order += 1

    # Store the last merge correctly
    last_unmerged = [seq for seq in new_distanceMatrix[0][1:] if seq not in flatten(list(merged.values()))]
    merged[order] = last_unmerged

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


def feng_doolittle_msa(sequences, score_matrix, gap_penalty, scoringType):
    """
    Implements the Feng-Doolittle multiple sequence alignment method.
    :param sequences: List of sequences to align
    :param score_matrix: Scoring matrix for pairwise alignments
    :param gap_penalty: Penalty for introducing gaps
    :return: Aligned sequences
    """
    # Step 1: Compute pairwise distance matrix using Needleman-Wunsch
    distance_dict, distance_matrix = buildDistanceMatrix(sequences, score_matrix, gap_penalty, scoringType)
    
    # Step 2: Construct a guide tree using hierarchical clustering (Fitch-Margoliash)
    _, merge_order = FitchMargoliash(distance_matrix, distance_dict)
    
    print(merge_order)

    # Step 3: Perform progressive alignment using the guide tree
    aligned_sequences = progressive_alignment(merge_order, sequences, gap_penalty, score_matrix)
    
    
    for row in aligned_sequences:
        print(row)
    
    return aligned_sequences[0], aligned_sequences[1]

def progressive_alignment(merge_order, sequences, gap_penalty, score_matrix):
    """
    Progressive alignment of sequences based on the guide tree.
    """
    aligned_sequences = copy.deepcopy(sequences)
    prev_aligned = {i: set() for i in range(len(sequences))}
    
    for order in merge_order:
        seq1, seq2 = merge_order[order]
        if isinstance(seq1, int): seq1 = [seq1]
        if isinstance(seq2, int): seq2 = [seq2]
        
        high_score = float('-inf')
        best_alignment = None
        
        for i in seq1:
            for j in seq2:
                aligned_pair, score = align_sequences(aligned_sequences[i-1], aligned_sequences[j-1], gap_penalty, score_matrix)
                if score > high_score:
                    high_score = score
                    best_alignment = (i, j, aligned_pair)
        
        print(i,j)
        i, j, (aligned_i, aligned_j) = best_alignment
        aligned_sequences[i-1] = aligned_i
        aligned_sequences[j-1] = aligned_j
        prev_aligned[i-1].add(j-1)
        prev_aligned[j-1].add(i-1)
        
        for k in prev_aligned[i-1]:
            if k != j-1:
                aligned_sequences[k] = add_gaps_to_match(aligned_sequences[i-1], aligned_sequences[k])
        for k in prev_aligned[j-1]:
            if k != i-1:
                aligned_sequences[k] = add_gaps_to_match(aligned_sequences[j-1], aligned_sequences[k])
        repl = [seq.replace('X', '-') for seq in aligned_sequences]
        aligned_sequences = repl

    return [seq.replace('X', '-') for seq in aligned_sequences], score

def align_sequences(seq1, seq2, gap_penalty, score_matrix):
    """
    Aligns two sequences using Needleman-Wunsch and returns the aligned sequences and alignment score.
    """
    aligned_seq, score = needleman_alg([seq1, seq2], gap_penalty, score_matrix)
    return aligned_seq, score

def add_gaps_to_match(ref_seq, seq):
    """
    Adds gaps to a sequence to match the reference sequence's gaps.
    """
    new_seq = []
    ref_index = 0
    for char in ref_seq:
        if char == '-':
            new_seq.append('-')
        else:
            new_seq.append(seq[ref_index])
            ref_index += 1
    return ''.join(new_seq)

def grabLabels(args):
    sequencelabels = []
    with open(args.input, "r") as file:
        genome = file.readlines()

    for row in genome:
        if row.startswith(">"):
            row = row.strip()
            sequencelabels.append(row)

    return sequencelabels

def flatten(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))  # Recursively flatten sublists
        else:
            result.append(item)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Feng-Doolittle Multiple Sequence Alignment")
    parser.add_argument("-i", "--input", required=True, type=str, help="File path to input.fna")
    parser.add_argument("-o", "--output", required=True, type=str, help="File path for output.fna")
    parser.add_argument("-s", "--score_matrix", required=True, type=str, help="File path for scoring matrix")
    args = parser.parse_args()
    
    sorted_sequences, genome = ReadAndStoreSequences(args)
    sorted_strings = outSortedSeqToFile(sorted_sequences, genome)
    score_matrix, gap_penalty = grabGapscore_matrix(args)
    scoringType = str(args.score_matrix)

    aligned_sequences, score = feng_doolittle_msa(sorted_strings, score_matrix, gap_penalty, scoringType)
    if(len(aligned_sequences[-1])<len(aligned_sequences[-2])):
        difference = len(aligned_sequences[-2])-len(aligned_sequences[-1])
        aligned_sequences[-1] += "-"*difference
    
    for row in aligned_sequences:
        print(row)

        
    #now writing to file for output:: 
    with open(args.output, "w") as output_file:
        for i in range(aligned_sequences):
            outputString = "%s; score=%d\n%s\n" % (sequenceLabels[i], msaScore, aligned_sequences[i])
            output_file.write(outputString)

