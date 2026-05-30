/**
 * Standalone interactive preview for the high-lift geometry engine.
 *
 * Vanilla TS, served by Vite at /highlift-preview.html. Deployment is a single
 * slider along the flap track (1-DOF); the remaining sliders shape the track and
 * the main cove. Phase B replaces this harness with the real /highlift route.
 */

import { DEFAULT_ESTOL_CONFIG } from './estolConfig';
import type { EstolConfig } from './estolConfig';
import { buildConfiguration, toPoints, trackPoint, trailingAxel } from './geometry';

const cfg: EstolConfig = structuredClone(DEFAULT_ESTOL_CONFIG);
let deploy = 0.3; // travel beyond stowed; sLead = axelSpacing + deploy

const VIEW = { xmin: -0.12, xmax: 1.55, ymin: -0.55, ymax: 0.38 };
const canvas = document.getElementById('view') as HTMLCanvasElement;
const ctx = canvas.getContext('2d')!;

function world(w: number, h: number) {
  const s = Math.min(w / (VIEW.xmax - VIEW.xmin), h / (VIEW.ymax - VIEW.ymin));
  const ox = (w - s * (VIEW.xmax - VIEW.xmin)) / 2;
  const oy = (h - s * (VIEW.ymax - VIEW.ymin)) / 2;
  return {
    x: (x: number) => ox + (x - VIEW.xmin) * s,
    y: (y: number) => h - (oy + (y - VIEW.ymin) * s),
  };
}

function fill(pts: { x: number; y: number }[], t: ReturnType<typeof world>, faceColor: string, edge: string) {
  ctx.beginPath();
  ctx.moveTo(t.x(pts[0].x), t.y(pts[0].y));
  for (const p of pts.slice(1)) ctx.lineTo(t.x(p.x), t.y(p.y));
  ctx.closePath();
  ctx.fillStyle = faceColor;
  ctx.fill();
  ctx.lineWidth = 1.4;
  ctx.strokeStyle = edge;
  ctx.stroke();
}

function render() {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = 540;
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const t = world(w, h);

  ctx.clearRect(0, 0, w, h);
  ctx.lineWidth = 1;
  ctx.strokeStyle = '#eef0f4';
  ctx.fillStyle = '#9aa0aa';
  ctx.font = '11px system-ui';
  for (let gx = 0; gx <= 1.4 + 1e-9; gx += 0.2) {
    ctx.beginPath(); ctx.moveTo(t.x(gx), t.y(VIEW.ymax)); ctx.lineTo(t.x(gx), t.y(VIEW.ymin)); ctx.stroke();
    ctx.fillText(gx.toFixed(1), t.x(gx) - 8, t.y(VIEW.ymin) - 4);
  }
  for (let gy = -0.4; gy <= 0.2 + 1e-9; gy += 0.1) {
    ctx.beginPath(); ctx.moveTo(t.x(VIEW.xmin), t.y(gy)); ctx.lineTo(t.x(VIEW.xmax), t.y(gy)); ctx.stroke();
    ctx.fillText(gy.toFixed(1), t.x(VIEW.xmin) + 2, t.y(gy) - 3);
  }

  // flap track (dashed), under the elements
  ctx.setLineDash([6, 5]);
  ctx.strokeStyle = '#9c5bd1';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(t.x(trackPoint(cfg.track, 0)[0]), t.y(trackPoint(cfg.track, 0)[1]));
  for (let i = 1; i <= 70; i++) {
    const p = trackPoint(cfg.track, (0.7 * i) / 70);
    ctx.lineTo(t.x(p[0]), t.y(p[1]));
  }
  ctx.stroke();
  ctx.setLineDash([]);

  const sLead = cfg.axelSpacing + deploy;
  const { main, vane, flap } = buildConfiguration(cfg, sLead);
  fill(toPoints(main), t, '#cfd8e6', '#1f3a68');
  fill(toPoints(vane), t, '#fbd4b4', '#9c3d1a');
  fill(toPoints(flap), t, '#cfe9c8', '#2c6b2c');

  // axels: the rigid bar between the two track carriages, on top
  const lead = trackPoint(cfg.track, sLead);
  const trail = trackPoint(cfg.track, trailingAxel(cfg.track, cfg.axelSpacing, sLead));
  ctx.strokeStyle = '#c0392b';
  ctx.lineWidth = 2.5;
  ctx.beginPath();
  ctx.moveTo(t.x(lead[0]), t.y(lead[1]));
  ctx.lineTo(t.x(trail[0]), t.y(trail[1]));
  ctx.stroke();
  ctx.fillStyle = '#c0392b';
  for (const a of [lead, trail]) {
    ctx.beginPath();
    ctx.arc(t.x(a[0]), t.y(a[1]), 4, 0, 2 * Math.PI);
    ctx.fill();
  }
}

function slider(label: string, min: number, max: number, step: number, get: () => number, set: (v: number) => void): HTMLElement {
  const row = document.createElement('div');
  row.className = 'row';
  const dp = step < 0.01 ? 3 : step < 1 ? 2 : 0;
  row.innerHTML = `<label>${label}</label>`;
  const input = document.createElement('input');
  input.type = 'range';
  input.min = String(min); input.max = String(max); input.step = String(step);
  input.value = String(get());
  const val = document.createElement('span');
  val.className = 'val';
  val.textContent = get().toFixed(dp);
  input.oninput = () => {
    set(parseFloat(input.value));
    val.textContent = parseFloat(input.value).toFixed(dp);
    render();
  };
  row.append(input, val);
  return row;
}

function naca(label: string, get: () => string, set: (v: string) => void): HTMLElement {
  const row = document.createElement('div');
  row.className = 'row';
  row.innerHTML = `<label>${label}</label>`;
  const input = document.createElement('input');
  input.type = 'text'; input.value = get(); input.maxLength = 4;
  input.oninput = () => /^\d{4}$/.test(input.value.trim()) && (set(input.value.trim()), render());
  row.append(input);
  return row;
}

function fieldset(legend: string, ...rows: HTMLElement[]): HTMLElement {
  const fs = document.createElement('fieldset');
  fs.innerHTML = `<legend>${legend}</legend>`;
  fs.append(...rows);
  return fs;
}

const { track: tr, mainCutouts: c } = cfg;
document.getElementById('controls')!.append(
  fieldset('Deployment',
    slider('deploy', 0, 0.45, 0.005, () => deploy, (v) => (deploy = v)),
  ),
  fieldset('Flap track',
    slider('start x', 0.4, 0.9, 0.005, () => tr.anchor[0], (v) => (tr.anchor[0] = v)),
    slider('start y', -0.2, 0.05, 0.005, () => tr.anchor[1], (v) => (tr.anchor[1] = v)),
    slider('axel gap', 0.05, 0.2, 0.005, () => cfg.axelSpacing, (v) => (cfg.axelSpacing = v)),
    slider('line angle°', -30, 10, 0.5, () => tr.angleDeg, (v) => (tr.angleDeg = v)),
    slider('linear len', 0.05, 0.3, 0.005, () => tr.linearLength, (v) => (tr.linearLength = v)),
    slider('arc radius', -0.25, -0.05, 0.005, () => tr.arcRadius, (v) => (tr.arcRadius = v)),
  ),
  fieldset('Main cove cutout',
    slider('lower cut x', 0.4, 0.62, 0.005, () => c.lowerCutX, (v) => (c.lowerCutX = v)),
    slider('upper lip x', 0.74, 0.95, 0.005, () => c.upperLipCutX, (v) => (c.upperLipCutX = v)),
    slider('cove vtx x', 0.54, 0.72, 0.005, () => c.coveVertexX, (v) => (c.coveVertexX = v)),
    slider('cove vtx y', 0.0, 0.1, 0.002, () => c.coveVertexY, (v) => (c.coveVertexY = v)),
    slider('fillet r', 0.05, 0.4, 0.005, () => c.coveFilletRadius, (v) => (c.coveFilletRadius = v)),
  ),
  fieldset('Elements (NACA 4-digit)',
    naca('vane', () => cfg.vane.naca, (v) => (cfg.vane.naca = v)),
    naca('aft flap', () => cfg.aftFlap.naca, (v) => (cfg.aftFlap.naca = v)),
  ),
);

render();
window.addEventListener('resize', render);
