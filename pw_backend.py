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
import os
import tempfile

app = Flask(__name__)
CORS(app)

@app.route('/', methods=['POST'])
def receive_sequences():
    # Check if file was uploaded
    if 'fastaFile' not in request.files:
        return jsonify({"error": "No file uploaded"}), 400
    
    file = request.files['fastaFile']
    sequenceType = request.form.get('sequenceType', 'DNA')
    
    if file.filename == '':
        return jsonify({"error": "No file selected"}), 400
    
    # Save uploaded file temporarily
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_file:
        file_content = file.read().decode('utf-8')
        temp_file.write(file_content)
        temp_file_path = temp_file.name
    
    try:
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

        sequenceList, labelToSeq, seqToLabel, titleList = readSequenceFile.main(temp_file_path)
        alignedSequences, msaScore = msa.main(sequenceList, gapPenalty, score_matrix)

        # Format the aligned sequences for better display
        formatted_sequences = "\n".join([f"Sequence {i+1}: {seq}" for i, seq in enumerate(alignedSequences)])
        
        return jsonify({
            "status": f"Sequences received and alignment completed successfully!\n\nAligned Sequences:\n{formatted_sequences}\n\nAlignment Score: {msaScore}"
        }), 200
    
    except Exception as e:
        return jsonify({"error": f"Error processing file, check sequence type"}), 500
    
    finally:
        # Clean up temporary file
        if os.path.exists(temp_file_path):
            os.unlink(temp_file_path)


if __name__ == '__main__':
    app.run(debug=True)
