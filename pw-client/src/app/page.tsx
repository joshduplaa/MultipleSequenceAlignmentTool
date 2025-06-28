'use client';

import { useState } from 'react';
import * as React from 'react';

export default function Home() {
  //defining object state (constructors) for backend
  const [fastaFile, setFastaFile] = useState<File | null>(null);
  const [sequenceType, setSequenceType] = useState<string>('DNA') //Store a selected sequence type (DNA or Protein)
  const [response, setResponse] = useState('');
  const [error, setError] = useState('');
  const [isLoading, setIsLoading] = useState(false);

  //event handler function for alignment type selection
  const handleSequenceType = (value: string) => {
    setSequenceType(value)
  }

  //event handler function for file upload
  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0];
    if (file) {
      setFastaFile(file);
    }
  }
  
  //Function to handle submit
  const handleSubmit = async () => {
    // Reset previous messages
    setResponse('');
    setError('');
    
    //Checks if the user has filled in the form correctly
    if (!sequenceType || !fastaFile) {
      setError('Make sure all options are selected and a FASTA file is uploaded!');
      return;
    }

    setIsLoading(true);

    try {
      // Create FormData to send file
      const formData = new FormData();
      formData.append('fastaFile', fastaFile);
      formData.append('sequenceType', sequenceType);

      //sends values to API, replace API route with secret when possible
      const res = await fetch('http://localhost:5000/', {
        method: 'POST',
        body: formData,
      });

      if (!res.ok) {
        // Handle HTTP errors (4xx, 5xx)
        const errorData = await res.json().catch(() => ({ error: 'Unknown error occurred' }));
        throw new Error(errorData.error || `Server error: ${res.status} ${res.statusText}`);
      }

      const data = await res.json();
      setResponse(data.status);
    } catch (err) {
      const errorMessage = err instanceof Error ? err.message : 'An unexpected error occurred';
      setError(errorMessage);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <>
      {/**Title */}
      <h1 className="title">Multiple Sequence Alignment Tool</h1>

      {/**Introduction */}
      <p className="intro">A web app for DNA and protein Multiple sequence alignment using Fitch-Margoliash (Distance Tree) and FengDoolittle (MSA) algorithms. Developed by Joshua Duplaa in Python for Dr. Rees&apos;t Bioinformatics course at TTU.</p>
      
      <div className='mainform'>

        {/**Sequence Type selection */}
        <div className="selector">
          <span>Input Sequence Type - </span>
          {['DNA', 'Protein'].map((option) => (
              <label key={option}>
                <input
                type="radio"
                name="sequence_type" //name of button group
                value={option}
                checked={sequenceType === option} //check if this option is alignType
                onChange={() => handleSequenceType(option)} //set the alignType value
              />
                <span>{option} </span>
              </label>
            ))}
        </div>
        <div className='input'>
            {/**Button to upload FASTA File */}
            <div className="file-upload">
              <label htmlFor="fasta-file" className="file-label">
                Choose FASTA File
              </label>
              <input
                id="fasta-file"
                type="file"
                accept=".fasta,.fna"
                onChange={handleFileUpload}
                className="file-input"
              />
              {fastaFile && (
                <p className="file-name">Selected file: {fastaFile.name}</p>
              )}
            </div>
        </div>
      </div>
      <button 
        className="bg-blue-500 text-white px-4 py-2 rounded disabled:bg-gray-400 disabled:cursor-not-allowed" 
        onClick={handleSubmit}
        disabled={isLoading}
      >
        {isLoading ? 'Processing...' : 'Submit'}
      </button>

      {/**Error message */}
      {error && (
        <div className="mt-4 p-4 bg-red-100 border border-red-400 text-red-700 rounded">
          <strong>Error:</strong> {error}
        </div>
      )}

      {/**Response from API */}
      {response && (
        <div className="mt-4 p-4 bg-green-100 border border-green-400 text-green-700 rounded">
          <strong>Success:</strong>
          <pre className="mt-2 whitespace-pre-wrap font-mono text-sm">{response}</pre>
        </div>
      )}
    </>
  );
}
