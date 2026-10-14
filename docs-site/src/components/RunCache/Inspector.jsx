import React, {useState} from 'react';

export default function RunCacheInspector() {
  const [current, setCurrent] = useState(true);
  const [converged, setConverged] = useState(true);
  const [trip, setTrip] = useState(1);
  const eligible = current && converged && trip === 1;
  return <div>
    <label><input type="checkbox" checked={current} onChange={e => setCurrent(e.target.checked)} /> Same known solver artifact</label>{' '}
    <label><input type="checkbox" checked={converged} onChange={e => setConverged(e.target.checked)} /> Successful converged result</label>{' '}
    <label>Requested upper trip x/c <input type="number" min="0" max="1" step="0.1" value={trip} onChange={e => setTrip(Number(e.target.value))} /></label>
    <p>{eligible ? 'Reuse the stored result.' : 'Compute and append a new run.'} The stored attempt remains in history either way. This example holds the other point conditions and settings fixed and uses a stored trip of 1.</p>
  </div>;
}
