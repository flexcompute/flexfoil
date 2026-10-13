import React, { useState } from 'react';
import { edgeState } from './math';

export default function CompressibilityInspector() {
  const [mach, setMach] = useState(.25);
  const [ui, setUi] = useState(1.2);
  const state = edgeState(ui, mach, .001, 1e6, 2.6);
  return <fieldset><legend>Compressible station</legend>
    <label style={{ display: 'block' }}>Freestream Mach: {mach.toFixed(2)} <input type="range" min="0" max="0.5" step="0.01" value={mach} onChange={e => setMach(Number(e.target.value))} /></label>
    <label style={{ display: 'block' }}>Incompressible speed / freestream: {ui.toFixed(2)} <input type="range" min="0.5" max="1.5" step="0.01" value={ui} onChange={e => setUi(Number(e.target.value))} /></label>
    <output aria-live="polite">{`Physical U = ${state.u.toFixed(6)}, local Mach = ${Math.sqrt(state.localMsq).toFixed(6)}, Hk = ${state.hk.toFixed(6)}, Rtheta = ${state.rtheta.toFixed(6)}`}</output>
    <p>Illustrative algebra only. The station uses theta = 0.001 chord, Re = 1,000,000 and H = 2.6. It does not run a coupled flow solution.</p>
  </fieldset>;
}
