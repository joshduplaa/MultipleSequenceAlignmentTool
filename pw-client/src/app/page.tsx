'use client';

import { useState } from 'react';
import * as React from 'react';

export default function Home() {
  //defining object state (constructors) for backend
  //const [fastaFile, setfastaFile] = useState(''); I don't think this linke is right
  const [sequenceType, setSequenceType] = useState<string>('DNA') //Store a selected sequence type (DNA or Protein)
  const [response, setResponse] = useState('');

  //event handler function for alignment type selection
  const handleSequenceType = (value: string) => {
    setSequenceType(value)
  }
  

  //Function to handle submit
  const handleSubmit = async () => {
    //Checks if the user has filled in the form correctly
    if ( !sequenceType || !fastaFile) {
      alert('Make sure all options are selected and the text boxes are not empty!');
      return;
    }
    //sends values to API, replace API route with secret when possible
    const res = await fetch('http://localhost:5000/', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ fastaFile, sequenceType}),
    });
    const data = await res.json();
    setResponse(data.status);
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

        </div>
      </div>
      <button className="bg-blue-500 text-white px-4 py-2 rounded" onClick={handleSubmit}>
        Submit
      </button>

      {/**Response from API */}
      {response && <p className="mt-4">{response}</p>}
    </>
  );
}
