import React, { useState } from 'react';
import coverage from '@site/src/generated/coverage.json';

export default function CoverageInspector() {
  const [query, setQuery] = useState('');
  const selected = coverage.filter(item => item.name.toLowerCase().includes(query.toLowerCase()));
  const count = selected.filter(item => item.page).length;
  return <fieldset><legend>Explore hand-written documentation coverage</legend>
    <label>Filter primitive names <input value={query} onChange={event => setQuery(event.target.value)} /></label>
    <p><output aria-live="polite">{count} of {selected.length} matching declarations have a hand-written primitive page.</output></p>
    <ul>{selected.slice(0, 15).map(item => <li key={`${item.path}:${item.line}`}>{item.name}: {item.page ? 'mapped' : 'needs a hand-written page'}</li>)}</ul>
    <p>The generated Coverage page lists every unmapped declaration. Rust source candidates may include declarations in private modules.</p>
  </fieldset>;
}
