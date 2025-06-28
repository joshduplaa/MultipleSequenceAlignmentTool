"""
=============================================================================
Title : pw_backend.py
Description : This is the backend flask app for my pairwise alignment app
Author : Joshua Duplaa
Date : 04/10/2025
Version : 1.0
Notes : this program has no requirements
Python Version: 3.11.3
=============================================================================
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
import util.msa as msa
import util.readSequenceFile as readSequenceFile

app = Flask(__name__)
CORS(app)

@app.route('/', methods=['POST'])
def receive_sequences():
    data = request.get_json()
    inputFile = data.get('inputFile')
    sequenceType = data.get('sequenceType', [])

    if sequenceType == "DNA":
        with open("nucleotide.mtx", "r") as file:
            score_matrix = file.readlines()
    elif sequenceType == "Protein":
        with open("BLOSUM50.mtx", "r") as file:
            score_matrix = file.readlines()

    scoreIndex = 0
    for row in score_matrix:
        row = row.strip("\n")
        score_matrix[scoreIndex] = row
        scoreIndex += 1
    score_matrix[0] = "0" + score_matrix[0]
    score_matrix = [row.split() for row in score_matrix]
    gapPenalty = int(score_matrix[3][len(score_matrix[3])-1])


    sequenceList, labelToSeq, seqToLabel, titleList = readSequenceFile.main(inputFile)
    alignedSequences, msaScore = msa.main(sequenceList, gapPenalty, score_matrix)

    return jsonify({"status": f"sequences received and alignment started \n {alignedSequences}, alignment score {msaScore}"}), 200


if __name__ == '__main__':
    app.run(debug=True)
