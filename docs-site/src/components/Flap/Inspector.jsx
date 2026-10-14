import React, {useState} from 'react';

export default function FlapInspector() {
  const [angle, setAngle] = useState(12);
  const [hinge, setHinge] = useState(.72);
  const radians = angle * Math.PI / 180;
  const x = hinge + (1 - hinge) * Math.cos(radians);
  const y = -(1 - hinge) * Math.sin(radians);
  return <div>
    <label>Hinge x/c <input type="number" min="0.01" max="0.99" step="0.01" value={hinge} onChange={event => setHinge(Number(event.target.value))} /></label>{' '}
    <label>Deflection (degrees) <input type="range" min="-30" max="30" value={angle} onChange={event => setAngle(Number(event.target.value))} />{angle}°</label>
    <svg viewBox="0 -0.25 1.05 0.5" style={{width:'100%',height:180}} role="img" aria-label="Rigid rotation of a trailing-edge point about the hinge">
      <line x1="0" y1="0" x2={hinge} y2="0" stroke="currentColor" strokeWidth="0.01" />
      <line x1={hinge} y1="0" x2={x} y2={-y} stroke="steelblue" strokeWidth="0.01" />
    </svg>
    <p>Trailing point: x/c={x.toFixed(5)}, y/c={y.toFixed(5)}. This shows the rotation step only; the real function also intersects and heals the airfoil contour.</p>
  </div>;
}
