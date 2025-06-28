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
import globalAlign
import localAlign

app = Flask(__name__)
CORS(app)

@app.route('/', methods=['POST'])
def receive_sequences():
    data = request.get_json()
    seq1 = data.get('seq1')
    seq2 = data.get('seq2')
    alignType = data.get('alignType', [])
    sequenceType = data.get('sequenceType', [])

    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print('User selected:', alignType)
    print('User selected:', sequenceType)

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

    if alignType == "Global":
        alignedSequences, alignmentScore = globalAlign.main([seq1, seq2], score_matrix)

    elif alignType == "Local":
        alignedSequences, alignmentScore = localAlign.main([seq1, seq2], score_matrix)

    return jsonify({"status": f"sequences received and alignment started \n {alignedSequences}, alignment score {alignmentScore}"}), 200


if __name__ == '__main__':
    app.run(debug=True)
