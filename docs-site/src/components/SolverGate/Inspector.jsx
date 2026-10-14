import React, { useState } from 'react';

export default function GateInspector() {
  const [reference, setReference] = useState(0.4432);
  const [actual, setActual] = useState(0.4433);
  const [tolerance, setTolerance] = useState(0.0002);
  const error = Math.abs(actual - reference);
  const passes = Number.isFinite(actual) && Number.isFinite(reference) && error <= tolerance;
  return <fieldset><legend>Drive the coefficient check</legend>
    {[['Reference CL', reference, setReference], ['Solver CL', actual, setActual], ['Absolute tolerance', tolerance, setTolerance]].map(([label, value, setter]) =>
      <label key={label} style={{ display: 'block', marginBlock: '0.5rem' }}>{label}{' '}
        <input type="number" step="0.00001" value={value} onChange={e => setter(e.target.valueAsNumber)} />
      </label>)}
    <output aria-live="polite">Absolute error: {error.toExponential(4)}. {passes ? 'Passes this coefficient check' : 'Fails this coefficient check'}.</output>
    <p>A complete case must also converge and pass every other coefficient, profile and binding check.</p>
  </fieldset>;
}
