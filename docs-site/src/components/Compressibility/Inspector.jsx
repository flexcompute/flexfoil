import React, { useState } from 'react';
import { edgeState, pressure, surfaceTransition } from './math';

export default function CompressibilityInspector({ mode = 'edge' }) {
  const [mach, setMach] = useState(.25);
  const [ui, setUi] = useState(1.2);
  const [itran, setItran] = useState(80);
  const state = edgeState(ui, mach, .001, 1e6, 2.6);
  if (mode === 'transition') return <fieldset><legend>Surface versus wake index</legend>
    <label>ITRAN (one based): {itran} <input type="range" min="1" max="84" value={itran} onChange={e => setItran(Number(e.target.value))} /></label>
    <p>Last surface index IBLTE (zero based): 80</p>
    <output aria-live="polite">{surfaceTransition(itran, 80) ? 'Report the interpolated surface location, including the final surface interval' : 'Transition lies in the wake; report 1.0 chord'}</output>
  </fieldset>;
  return <fieldset><legend>{mode === 'pressure' ? 'Karman–Tsien pressure' : 'Compressible station'}</legend>
    <label style={{ display: 'block' }}>Freestream Mach: {mach.toFixed(2)} <input type="range" min="0" max="0.5" step="0.01" value={mach} onChange={e => setMach(Number(e.target.value))} /></label>
    <label style={{ display: 'block' }}>Incompressible speed / freestream: {ui.toFixed(2)} <input type="range" min="0.5" max="1.5" step="0.01" value={ui} onChange={e => setUi(Number(e.target.value))} /></label>
    <output aria-live="polite">{mode === 'pressure' ? `Cp = ${pressure(ui, mach).toFixed(6)} (incompressible Cp = ${(1-ui*ui).toFixed(6)})` : `Physical U = ${state.u.toFixed(6)}, local Mach = ${Math.sqrt(state.localMsq).toFixed(6)}, Hk = ${state.hk.toFixed(6)}, Rtheta = ${state.rtheta.toFixed(6)}`}</output>
    <p>Illustrative algebra only. The station uses theta = 0.001 chord, Re = 1,000,000 and H = 2.6. It does not run a coupled flow solution.</p>
  </fieldset>;
}
