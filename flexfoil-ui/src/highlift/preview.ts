/**
 * Standalone interactive high-lift design tool (geometry + per-element aero).
 *
 * Vanilla TS, served by Vite at /highlift-preview.html. Deployment is a single
 * slider along the flap track (1-DOF); the other sliders shape the track, cove,
 * hinge and NACA elements. The drag polar and Cp plots come from the per-element
 * independent solver and re-run automatically as inputs change. The polar run is
 * async + cancellable: it greys out while computing and restarts if geometry or
 * polar parameters change mid-run.
 */

// Prebuilt browser bundle (self-contained) — the plotly.js source entry pulls
// in node polyfills (buffer/) that break Vite's dep pre-bundle.
import Plotly from 'plotly.js/dist/plotly-basic.min.js';
import { DEFAULT_HIGH_LIFT_AIRFOIL } from './estolConfig';
import type { HighLiftAirfoil } from './estolConfig';
import { buildConfiguration, toPoints, trackPoint, trailingAxel } from './geometry';
// The solver (and the WASM layer it imports) is loaded lazily so the geometry
// view still works when the WASM package isn't built. Types are erased.
import type { ElementName } from './solve';

const cfg: HighLiftAirfoil = structuredClone(DEFAULT_HIGH_LIFT_AIRFOIL);
let alphaDeg = 4;
const COLORS: Record<ElementName, string> = { main: '#1f3a68', vane: '#9c3d1a', flaperon: '#2c6b2c' };

const VIEW = { xmin: -0.1, xmax: 1.2, ymin: -0.3, ymax: 0.2 };
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
  const h = 360;
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const t = world(w, h);

  ctx.clearRect(0, 0, w, h);
  ctx.lineWidth = 1;
  ctx.strokeStyle = '#eef0f4';
  ctx.fillStyle = '#9aa0aa';
  ctx.font = '11px system-ui';
  for (let gx = 0; gx <= 1.2 + 1e-9; gx += 0.2) {
    ctx.beginPath(); ctx.moveTo(t.x(gx), t.y(VIEW.ymax)); ctx.lineTo(t.x(gx), t.y(VIEW.ymin)); ctx.stroke();
    ctx.fillText(gx.toFixed(1), t.x(gx) - 8, t.y(VIEW.ymin) - 4);
  }
  for (let gy = -0.3; gy <= 0.2 + 1e-9; gy += 0.1) {
    ctx.beginPath(); ctx.moveTo(t.x(VIEW.xmin), t.y(gy)); ctx.lineTo(t.x(VIEW.xmax), t.y(gy)); ctx.stroke();
    ctx.fillText(gy.toFixed(1), t.x(VIEW.xmin) + 2, t.y(gy) - 3);
  }

  // flap track (dashed), under the elements
  ctx.setLineDash([6, 5]);
  ctx.strokeStyle = '#9c5bd1';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(t.x(trackPoint(cfg.design.track, 0)[0]), t.y(trackPoint(cfg.design.track, 0)[1]));
  for (let i = 1; i <= 70; i++) {
    const p = trackPoint(cfg.design.track, (0.7 * i) / 70);
    ctx.lineTo(t.x(p[0]), t.y(p[1]));
  }
  ctx.stroke();
  ctx.setLineDash([]);

  const sLead = cfg.design.axelSpacing + cfg.operation.deploy;
  const { main, vane, flaperon } = buildConfiguration(cfg);
  fill(toPoints(main), t, '#cfd8e6', '#1f3a68');
  fill(toPoints(vane), t, '#fbd4b4', '#9c3d1a');
  fill(toPoints(flaperon), t, '#cfe9c8', '#2c6b2c');

  // axels: the rigid bar between the two track carriages, on top
  const lead = trackPoint(cfg.design.track, sLead);
  const trail = trackPoint(cfg.design.track, trailingAxel(cfg.design.track, cfg.design.axelSpacing, sLead));
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

// --- control builders ---
function slider(
  label: string, min: number, max: number, step: number,
  get: () => number, set: (v: number) => void, after: () => void = refresh,
): HTMLElement {
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
    after();
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
  input.oninput = () => /^\d{4}$/.test(input.value.trim()) && (set(input.value.trim()), refresh());
  row.append(input);
  return row;
}

function fieldset(legend: string, ...rows: HTMLElement[]): HTMLElement {
  const fs = document.createElement('fieldset');
  fs.innerHTML = `<legend>${legend}</legend>`;
  fs.append(...rows);
  return fs;
}

/** A read-only row (e.g. the fixed main airfoil at level 1). */
function info(label: string, value: string): HTMLElement {
  const row = document.createElement('div');
  row.className = 'row';
  row.innerHTML = `<label>${label}</label><span class="val" style="flex:1;text-align:left;color:#c8c8d4">${value}</span>`;
  return row;
}

// --- aero: solving runs in a Web Worker (off the main thread) ---
const worker = new Worker(new URL('./solveWorker.ts', import.meta.url), { type: 'module' });
const setPolarBusy = (b: boolean) => document.getElementById('polar')!.classList.toggle('busy', b);
let cpReq = 0;
let polarReq = 0;

worker.onmessage = (e: MessageEvent) => {
  const m = e.data;
  if (m.kind === 'cp' && m.id === cpReq) drawCp(m.res);
  else if (m.kind === 'polar' && m.id === polarReq) { drawPolar(m.rows); setPolarBusy(false); }
};
worker.onerror = () => {
  for (const id of ['polar', 'cp']) {
    document.getElementById(id)!.innerHTML =
      '<div style="padding:16px;color:#444;font:13px system-ui">Aero solver unavailable — build the WASM package (see README).</div>';
  }
};

function drawCp(res: { name: ElementName; cpX: number[]; cp: number[] }[]) {
  const traces = res.map((r) => ({ x: r.cpX, y: r.cp, name: r.name, mode: 'lines', line: { color: COLORS[r.name] } }));
  Plotly.react('cp', traces, {
    title: { text: `Cp at α = ${alphaDeg.toFixed(1)}°`, font: { size: 13 } },
    xaxis: { title: 'x/c' },
    yaxis: { title: 'Cp', autorange: 'reversed' },
    margin: { t: 30, r: 10, b: 40, l: 48 },
    legend: { orientation: 'h', y: -0.2 },
  }, { displayModeBar: false, responsive: true });
}

function drawPolar(rows: Record<ElementName, { alpha: number; cl: number }[]>) {
  const traces = (['main', 'vane', 'flaperon'] as ElementName[]).map((name) => ({
    x: rows[name].map((p) => p.alpha),
    y: rows[name].map((p) => p.cl),
    name, mode: 'lines+markers', line: { color: COLORS[name] }, marker: { size: 4 },
  }));
  Plotly.react('polar', traces, {
    title: { text: 'Cl vs α (coupled inviscid)', font: { size: 13 } },
    xaxis: { title: 'α (deg)' },
    yaxis: { title: 'Cl' },
    margin: { t: 30, r: 10, b: 40, l: 48 },
    legend: { orientation: 'h', y: -0.2 },
  }, { displayModeBar: false, responsive: true });
}

function requestCp() {
  worker.postMessage({ id: ++cpReq, kind: 'cp', cfg, alphaDeg });
}
function requestPolar() {
  setPolarBusy(true);
  worker.postMessage({ id: ++polarReq, kind: 'polar', cfg });
}

let cpTimer = 0;
let polarTimer = 0;
function scheduleCp() { clearTimeout(cpTimer); cpTimer = window.setTimeout(requestCp, 150); }
function schedulePolar() { clearTimeout(polarTimer); polarTimer = window.setTimeout(requestPolar, 200); }

function refresh() {
  render();
  scheduleCp();
  schedulePolar();
}

// --- geometry controls, grouped by the 3-level hierarchy ---
const m = cfg.main;          // L1
const d = cfg.design;        // L2
const op = cfg.operation;    // L3
const tr = d.track;
const c = d.cove;
const pivot = d.flaperonHinge.pivot;
document.getElementById('controls')!.append(
  // ── Level 1: main airfoil (foundational; fixed for now) ──
  fieldset('① Main airfoil',
    info('airfoil', 'LS(1)-0417 (fixed)'),
    slider('blunt TE', 0.0, 0.02, 0.001, () => m.bluntThickness, (v) => (m.bluntThickness = v)),
  ),
  // ── Level 2: high-lift element design ──
  fieldset('② High-lift design',
    fieldset('Main cove cutout',
      slider('lower cut x', 0.4, 0.62, 0.005, () => c.lowerCutX, (v) => (c.lowerCutX = v)),
      slider('upper lip x', 0.74, 0.95, 0.005, () => c.upperLipCutX, (v) => (c.upperLipCutX = v)),
      slider('cove vtx x', 0.54, 0.72, 0.005, () => c.coveVertexX, (v) => (c.coveVertexX = v)),
      slider('cove vtx y', 0.0, 0.1, 0.002, () => c.coveVertexY, (v) => (c.coveVertexY = v)),
      slider('fillet r', 0.05, 0.4, 0.005, () => c.coveFilletRadius, (v) => (c.coveFilletRadius = v)),
    ),
    fieldset('Vane & flaperon (NACA 4-digit)',
      naca('vane', () => d.vane.naca, (v) => (d.vane.naca = v)),
      naca('flaperon', () => d.flaperon.naca, (v) => (d.flaperon.naca = v)),
    ),
    fieldset('Flap track',
      slider('start x', 0.4, 0.9, 0.005, () => tr.anchor[0], (v) => (tr.anchor[0] = v)),
      slider('start y', -0.2, 0.05, 0.005, () => tr.anchor[1], (v) => (tr.anchor[1] = v)),
      slider('axel gap', 0.05, 0.2, 0.005, () => d.axelSpacing, (v) => (d.axelSpacing = v)),
      slider('line angle°', -30, 10, 0.5, () => tr.angleDeg, (v) => (tr.angleDeg = v)),
      slider('linear len', 0.05, 0.3, 0.005, () => tr.linearLength, (v) => (tr.linearLength = v)),
      slider('arc radius', -0.25, -0.05, 0.005, () => tr.arcRadius, (v) => (tr.arcRadius = v)),
    ),
    fieldset('Flaperon hinge (pivot)',
      slider('pivot x', 0.6, 1.05, 0.005, () => pivot[0], (v) => (pivot[0] = v)),
      slider('pivot y', -0.15, 0.1, 0.005, () => pivot[1], (v) => (pivot[1] = v)),
    ),
  ),
  // ── Level 3: deployment & flaperon angle (operating point) ──
  fieldset('③ Deployment',
    slider('deploy', 0, 0.45, 0.005, () => op.deploy, (v) => (op.deploy = v)),
    slider('flaperon angle°', -60, 30, 0.5, () => op.flaperonAngleDeg, (v) => (op.flaperonAngleDeg = v)),
  ),
);

// Polar window has no controls in the inviscid phase (no Re/Ncrit/drag yet);
// Re/Ncrit and a true drag polar return with the viscous phase.

// --- Cp window control (α only updates Cp, not the polar) ---
document.getElementById('cp-controls')!.append(
  slider('α (deg)', -5, 15, 0.5, () => alphaDeg, (v) => (alphaDeg = v), scheduleCp),
);

render();
window.addEventListener('resize', render);
requestCp();
requestPolar();
