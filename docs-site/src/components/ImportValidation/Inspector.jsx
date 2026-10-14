import React, {useState} from 'react';

export default function ImportValidationInspector() {
  const [line, setLine] = useState('0 Infinity');
  const parts = line.trim().replace(/,/g,' ').split(/\s+/);
  const x = Number(parts[0]);
  const y = Number(parts[1]);
  const numericLike = /^[+-]?(?:\d|\.\d|NaN\b|Inf(?:inity)?\b)/i.test(parts[0]);
  const definitePair = numericLike && (!Number.isNaN(y) || /^[+-]?NaN$/i.test(parts[1]));
  const result = Number.isFinite(x) && Number.isFinite(y)
    ? 'Both coordinate values are finite; the full parser still checks scale, duplicates and group format.'
    : definitePair || /^[+-]?(?:NaN|Inf(?:inity)?)$/i.test(parts[0])
      ? 'Reject this non-finite/malformed coordinate row, even before a header.'
      : 'This may be a text/numeric airfoil name when a header is allowed.';
  return <fieldset><legend>First coordinate pair or name</legend>
    <label>Line <input value={line} onChange={event=>setLine(event.target.value)} /></label>
    <p>{result}</p><p>Illustrates the finite-coordinate decision only; this does not run a complete DAT import.</p>
  </fieldset>;
}
