/**
 * Standalone interactive high-lift design tool (geometry + RANS flow field).
 *
 * Vanilla TS, served by Vite at /highlift-preview.html. The operation is a deploy
 * slider along the flap track (1-DOF) + a flap angle; the other sliders shape the
 * track, cove, vane and flaperon. Geometry drawing is live. Aero comes from the RANS
 * bridge (rans_server.py): "Run RANS" sweeps α and plots a clickable CD–CL polar (one
 * per operation, overlaid); the selected point's center-span flow field is drawn as a
 * two-pass LIC overlay (see lic.ts) behind the airfoil canvas.
 */

// Prebuilt browser bundle (self-contained) — the plotly.js source entry pulls
// in node polyfills (buffer/) that break Vite's dep pre-bundle.
import Plotly from 'plotly.js/dist/plotly-basic.min.js';
import { DEFAULT_HIGH_LIFT_AIRFOIL } from './estolConfig';
import type { HighLiftAirfoil } from './estolConfig';
import { buildConfiguration, toPoints, trackPoint, trailingAxel } from './geometry';
import { LicView } from './lic';
import type { FlowMesh } from './lic';

// This is a stateful vanilla entry (module-level state + appended DOM + Plotly
// listeners), so hot-swapping would stack duplicate handlers and run twice. Force a
// full reload on any change instead.
if ((import.meta as any).hot) (import.meta as any).hot.accept(() => location.reload());

const cfg: HighLiftAirfoil = structuredClone(DEFAULT_HIGH_LIFT_AIRFOIL);

const VIEW = { xmin: -0.1, xmax: 1.2, ymin: -0.4, ymax: 0.2 };
const canvas = document.getElementById('view') as HTMLCanvasElement;
const ctx = canvas.getContext('2d')!;
const VIEW_H = 360;

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
  const h = VIEW_H;
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
  for (let gy = VIEW.ymin; gy <= VIEW.ymax + 1e-9; gy += 0.1) {
    ctx.beginPath(); ctx.moveTo(t.x(VIEW.xmin), t.y(gy)); ctx.lineTo(t.x(VIEW.xmax), t.y(gy)); ctx.stroke();
    ctx.fillText(gy.toFixed(1), t.x(VIEW.xmin) + 2, t.y(gy) - 3);
  }

  // flap track (dashed) + axels — only while the Flap-track pane is open.
  if (trackPane.open) {
    ctx.setLineDash([6, 5]);
    ctx.strokeStyle = '#9c5bd1';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    const trackLen = cfg.design.track.linearLength + cfg.design.axelSpacing + cfg.design.track.length;
    ctx.moveTo(t.x(trackPoint(cfg.design.track, 0)[0]), t.y(trackPoint(cfg.design.track, 0)[1]));
    for (let i = 1; i <= 70; i++) {
      const p = trackPoint(cfg.design.track, (trackLen * i) / 70);
      ctx.lineTo(t.x(p[0]), t.y(p[1]));
    }
    ctx.stroke();
    ctx.setLineDash([]);
    // axels: the rigid bar between the two track carriages
    const sLead = cfg.design.axelSpacing + cfg.operation.deploy;
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

  const { main, vane, vaneControls, flaperon, flaperonControls, flaperonPivot } = buildConfiguration(cfg);
  fill(toPoints(main), t, '#cfd8e6', '#1f3a68');
  fill(toPoints(vane), t, '#fbd4b4', '#9c3d1a');
  fill(toPoints(flaperon), t, '#cfe9c8', '#2c6b2c');

  // B-spline control points / polygon, shown only while the matching pane is open.
  const drawControlNet = (pts: number[][], color: string, closed: boolean) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(t.x(pts[0][0]), t.y(pts[0][1]));
    for (const p of pts.slice(1)) ctx.lineTo(t.x(p[0]), t.y(p[1]));
    if (closed) ctx.closePath();
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = color;
    ctx.font = '10px system-ui';
    pts.forEach((p, i) => {
      ctx.beginPath();
      ctx.arc(t.x(p[0]), t.y(p[1]), 3.5, 0, 2 * Math.PI);
      ctx.fill();
      ctx.fillText(`P${i}`, t.x(p[0]) + 5, t.y(p[1]) - 5);
    });
  };
  if (vanePane.open) drawControlNet(vaneControls, '#1b8a8a', true);        // closed loop (teal)
  if (flaperonPane.open) drawControlNet(flaperonControls, '#7a4fc0', false); // open nose (purple)

  // flaperon hinge pivot — only while the Flaperon-hinge pane is open.
  if (hingePane.open) {
    const [px, py] = [t.x(flaperonPivot[0]), t.y(flaperonPivot[1])];
    ctx.strokeStyle = '#e08a1e';
    ctx.fillStyle = '#e08a1e';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(px, py, 5, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(px - 8, py); ctx.lineTo(px + 8, py);
    ctx.moveTo(px, py - 8); ctx.lineTo(px, py + 8);
    ctx.stroke();
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

/** A collapsible section (collapsed by default). Nests for the design sub-boxes. */
function section(title: string, ...rows: HTMLElement[]): HTMLDetailsElement {
  const d = document.createElement('details');           // no `open` attr ⇒ collapsed
  d.innerHTML = `<summary>${title}</summary>`;
  d.append(...rows);
  return d;
}

/** A read-only row (e.g. the fixed main airfoil at level 1). */
function info(label: string, value: string): HTMLElement {
  const row = document.createElement('div');
  row.className = 'row';
  row.innerHTML = `<label>${label}</label><span class="val" style="flex:1;text-align:left;color:#c8c8d4">${value}</span>`;
  return row;
}

// --- RANS flow field (LIC overlay) + clickable CD–CL polar, via the local bridge ---
// The GPU solver can't run in the browser; a local `rans_server.py` (user's shell)
// runs the fast pipeline and returns forces + a center-span slice mesh per α. The
// LIC layer draws that flow field behind the 2D airfoil canvas; the CD–CL polar is
// clickable to switch which point's flow field is shown.
const licCanvas = document.getElementById('view-lic') as HTMLCanvasElement;
const licView = new LicView(licCanvas);

interface PolarPoint { alpha: number; CL: number; CD: number; flow: FlowMesh | null; }
// One polar per operation point (deploy + flap angle), swept over α. They overlay so
// different deployments can be compared; the design geometry is shared across them.
interface Polar { deploy: number; flapAngle: number; quality: string; color: string; points: PolarPoint[]; }
// `polars` holds the sweeps; `current` is the cursor into them — the one point shown in
// the LIC and ring-highlighted on the polar. These two are the entire view state: set
// `current` (or clear polars) and call renderFlow()/drawPolar(); everything follows.
let polars: Polar[] = [];
let current: { polar: number; point: number } | null = null;
let ransBusy = false;
let meshQuality: 'fast' | 'accurate' = 'accurate';   // mesh refinement for the next sweep

const ALPHA_SWEEP = [-2, 0, 2, 4, 6, 8];
const POLAR_COLORS = ['#1f3a68', '#9c3d1a', '#2c6b2c', '#7a4fc0', '#b8860b', '#1b8a8a'];

function setStatus(msg: string): void {
  document.getElementById('rans-status')!.textContent = msg;
}

/** Size the LIC canvas to the airfoil canvas and align its camera to the 2D view. */
function syncLicCamera(): void {
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = VIEW_H;
  licCanvas.style.width = w + 'px';
  licCanvas.style.height = h + 'px';
  licView.resize(Math.round(w * dpr), Math.round(h * dpr));
  // world rect covering the full pixel buffer (mirrors world()'s letterbox math)
  const s = Math.min(w / (VIEW.xmax - VIEW.xmin), h / (VIEW.ymax - VIEW.ymin));
  const ox = (w - s * (VIEW.xmax - VIEW.xmin)) / 2;
  const oy = (h - s * (VIEW.ymax - VIEW.ymin)) / 2;
  licView.setCamera(VIEW.xmin - ox / s, VIEW.xmax + ox / s, VIEW.ymin - oy / s, VIEW.ymax + oy / s);
}

/** Draw the current point's flow field on the LIC overlay (or clear it if no point is
 *  selected). `current` alone decides what's shown. */
function renderFlow(): void {
  syncLicCamera();
  const flow = current && polars[current.polar]?.points[current.point]?.flow;
  if (flow) { licView.setFlow(flow); licView.render(); }
  else licView.clear();
}

function drawPolar(): void {
  const gd = document.getElementById('polar')!;
  const base = {
    xaxis: { title: 'CD' }, yaxis: { title: 'CL' }, margin: { t: 30, r: 10, b: 56, l: 52 },
    legend: { orientation: 'h', y: -0.2, font: { size: 11 } },
  };
  if (!polars.length) {
    Plotly.react('polar', [], { title: { text: 'CD vs CL (RANS) — click “Run RANS”', font: { size: 13 } }, ...base },
      { displayModeBar: false, responsive: true });
    return;
  }
  const traces: any[] = polars.map((pl) => ({
    x: pl.points.map((p) => p.CD), y: pl.points.map((p) => p.CL),
    name: `d=${pl.deploy.toFixed(2)} δ=${pl.flapAngle.toFixed(0)}° ${pl.quality[0]}`,
    text: pl.points.map((p) => `α=${p.alpha}°`), mode: 'lines+markers',
    line: { color: pl.color }, marker: { size: 8, color: pl.color },
    hovertemplate: '%{fullData.name}<br>%{text}<br>CD=%{x:.4f}  CL=%{y:.3f}<extra></extra>',
  }));
  const cur = current && polars[current.polar]?.points[current.point];
  if (cur) traces.push({                                  // ring on the point shown in the LIC
    x: [cur.CD], y: [cur.CL], mode: 'markers', showlegend: false, hoverinfo: 'skip',
    marker: { size: 16, symbol: 'circle-open', color: '#111', line: { width: 3 } },
  });
  Plotly.react('polar', traces, { title: { text: 'CD vs CL (RANS) — click a point', font: { size: 13 } }, ...base },
    { displayModeBar: false, responsive: true });
  bindPolarClick(gd);
}

let polarClickBound = false;
/** Bind the polar's click handler once (it reads live module state, so it needn't rebind). */
function bindPolarClick(gd: HTMLElement): void {
  if (polarClickBound) return;
  polarClickBound = true;
  (gd as any).on('plotly_click', (ev: any) => {
    const pt = ev.points[0];
    if (pt.curveNumber >= polars.length) return;          // the highlight marker
    selectPoint(pt.curveNumber, pt.pointNumber ?? pt.pointIndex ?? 0);
  });
}

/** Click a polar point: show its flow field and move the operation (deploy + flap
 *  angle) to that polar's state, so the airfoil view matches what's shown. */
function selectPoint(pi: number, ki: number): void {
  if (!polars[pi]?.points[ki]) return;
  current = { polar: pi, point: ki };
  setOperation(polars[pi].deploy, polars[pi].flapAngle);  // moves the flap geometry + the sliders
  renderFlow();
  drawPolar();                                            // reposition the highlight
}

/** The deployed configuration as the rans bridge's element list ([[x,y]…] contours). */
function ransElements() {
  const { main, vane, flaperon } = buildConfiguration(cfg);
  return [
    { name: 'main', contour: main },
    { name: 'vane', contour: vane },
    { name: 'flaperon', contour: flaperon },
  ];
}

const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

async function runRans(): Promise<void> {
  if (ransBusy) return;
  ransBusy = true;
  const btn = document.getElementById('run-rans') as HTMLButtonElement;
  btn.disabled = true;
  document.getElementById('polar')!.classList.add('busy');
  const elements = ransElements();
  // A polar is one operation point (current deploy + flap angle) swept over α. Replace
  // an existing polar at the same operation, else add a new (differently-coloured) one.
  const deploy = op.deploy;
  const flapAngle = op.flaperonAngleDeg;
  const quality = meshQuality;
  const key = (pl: { deploy: number; flapAngle: number; quality: string }) =>
    `${pl.deploy.toFixed(3)}|${pl.flapAngle.toFixed(1)}|${pl.quality}`;
  let idx = polars.findIndex((pl) => key(pl) === key({ deploy, flapAngle, quality }));
  const color = idx >= 0 ? polars[idx].color : POLAR_COLORS[polars.length % POLAR_COLORS.length];
  const polar: Polar = { deploy, flapAngle, quality, color, points: [] };
  if (idx >= 0) polars[idx] = polar; else { idx = polars.length; polars.push(polar); }
  current = null;
  drawPolar();
  // One unsteady-as-steady run marches through every α (warm-started); points stream
  // back as each physical step converges, so the polar fills in live.
  setStatus(`Starting RANS α sweep (${ALPHA_SWEEP.length} points, one run)…`);
  try {
    const start = await fetch('/api/rans/sweep', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ elements, alphas: ALPHA_SWEEP, quality }),
    }).then((r) => r.json());
    if (start.error) throw new Error(start.error);
    const { jobId, n } = start;
    let fails = 0;
    let shown = 0;
    for (;;) {
      await sleep(1500);
      let st: any;
      try {
        st = await fetch(`/api/rans/sweep/status?job=${jobId}`).then((r) => r.json());
      } catch (e) {
        if (++fails > 5) throw e;                 // tolerate transient poll failures
        continue;
      }
      fails = 0;
      polar.points = (st.points ?? []).map((p: any) => ({
        alpha: p.alpha, CL: p.CL, CD: p.CD, flow: p.flowField ?? null,
      }));
      if (polar.points.length > shown) {          // a new α landed → preview + highlight it
        shown = polar.points.length;
        current = { polar: idx, point: shown - 1 };
        renderFlow();
      }
      drawPolar();
      if (st.done) {
        setStatus(st.error
          ? `${st.error} — ${polar.points.length} converged point(s) shown.`
          : `Done — ${polar.points.length} points. Click any polar point to view its flow field.`);
        break;
      }
      const stage = polar.points.length ? `${polar.points.length}/${n} α` : (st.stage ?? 'meshing…');
      setStatus(`Running RANS sweep — ${stage}`);
    }
  } catch (err) {
    setStatus(`RANS error: ${(err as Error).message}. Is rans_server.py running in your shell?`);
  } finally {
    document.getElementById('polar')!.classList.remove('busy');
    btn.disabled = false;
    ransBusy = false;
  }
}

/** A design (geometry) change invalidates every polar — they were swept on the old
 *  shape. An operation-only change (deploy / flap angle) keeps the polars (the user can
 *  run another to overlay it). Either way the shown flow + highlight are now stale. */
function refresh() {
  render();
  if (ransBusy) return;
  polars = [];
  invalidateView();
}

function refreshOperation() {
  render();
  if (ransBusy) return;
  invalidateView();
}

/** Drop the displayed flow + highlight and repaint. */
function invalidateView() {
  current = null;
  renderFlow();
  drawPolar();
}

/** Move the operation to a state and reflect it in the two operation sliders. */
function setOperation(deploy: number, flapAngle: number): void {
  op.deploy = deploy;
  op.flaperonAngleDeg = flapAngle;
  syncSlider(deployRow, deploy, 3);
  syncSlider(flapRow, flapAngle, 1);
  render();
}

function syncSlider(row: HTMLElement, value: number, dp: number): void {
  (row.querySelector('input') as HTMLInputElement).value = String(value);
  (row.querySelector('.val') as HTMLElement).textContent = value.toFixed(dp);
}

// --- geometry controls, grouped by the 3-level hierarchy ---
const m = cfg.main;          // L1
const d = cfg.design;        // L2
const op = cfg.operation;    // L3
const tr = d.track;
const c = d.cove;
const f = d.flaperon;
const pivot = d.flaperonHinge.pivot;

// The flaperon sub-pane is held in a variable so render() can show the B-spline
// control points only while it is open (and re-render on toggle).
const flaperonPane = section('Flaperon shape',
  slider('upper cut x', 0.5, 0.95, 0.005, () => f.upperCutX, (v) => (f.upperCutX = v)),
  slider('lower cut x', 0.5, 0.95, 0.005, () => f.lowerCutX, (v) => (f.lowerCutX = v)),
  slider('nose tan up', 0.0, 0.15, 0.005, () => f.noseTangentUpper, (v) => (f.noseTangentUpper = v)),
  slider('nose tan low', 0.0, 0.15, 0.005, () => f.noseTangentLower, (v) => (f.noseTangentLower = v)),
  slider('nose tip x', 0.4, 0.9, 0.005, () => f.noseTip[0], (v) => (f.noseTip[0] = v)),
  slider('nose tip y', -0.1, 0.15, 0.005, () => f.noseTip[1], (v) => (f.noseTip[1] = v)),
);
flaperonPane.addEventListener('toggle', render);   // show/hide control points immediately

// Vane shape pane (its B-spline control points show while it is open, like the flaperon).
const v = d.vane;
const vanePane = section('Vane shape',
  slider('te x', 0.55, 0.85, 0.005, () => v.te[0], (x) => (v.te[0] = x)),
  slider('te y', -0.05, 0.12, 0.005, () => v.te[1], (x) => (v.te[1] = x)),
  ...v.points.flatMap((p, i) => [
    slider(`p${i + 1} x`, 0.45, 0.85, 0.005, () => p[0], (x) => (p[0] = x)),
    slider(`p${i + 1} y`, -0.12, 0.12, 0.005, () => p[1], (x) => (p[1] = x)),
  ]),
);
vanePane.addEventListener('toggle', render);

// Operation sliders (deploy + flap angle) use refreshOperation: changing them keeps
// the existing polars (so a new one can overlay) rather than clearing like a design edit.
// Deploy's max is the track `length`, kept in sync when the length slider moves.
const deployRow = slider('deploy', 0, tr.length, 0.005, () => op.deploy, (val) => (op.deploy = val), refreshOperation);
const flapRow = slider('flaperon angle°', -60, 30, 0.5, () => op.flaperonAngleDeg, (v) => (op.flaperonAngleDeg = v), refreshOperation);
const deployInput = deployRow.querySelector('input') as HTMLInputElement;
const deployVal = deployRow.querySelector('.val') as HTMLElement;
function setDeployMax(maxv: number): void {
  deployInput.max = String(maxv);
  if (op.deploy > maxv) {
    op.deploy = maxv;
    deployInput.value = String(maxv);
    deployVal.textContent = maxv.toFixed(3);
  }
}

// Track + hinge panes are held in variables so render() can show the track/axels
// and the flaperon pivot only while the respective pane is open.
const trackPane = section('Flap track',
  slider('start x', 0.4, 0.9, 0.005, () => tr.anchor[0], (v) => (tr.anchor[0] = v)),
  slider('start y', -0.2, 0.05, 0.005, () => tr.anchor[1], (v) => (tr.anchor[1] = v)),
  slider('axel gap', 0.05, 0.2, 0.005, () => d.axelSpacing, (v) => (d.axelSpacing = v)),
  slider('line angle°', -30, 10, 0.5, () => tr.angleDeg, (v) => (tr.angleDeg = v)),
  slider('linear len', 0.05, 0.3, 0.005, () => tr.linearLength, (v) => (tr.linearLength = v)),
  slider('arc radius', -0.25, -0.05, 0.005, () => tr.arcRadius, (v) => (tr.arcRadius = v)),
  slider('length', 0.05, 0.5, 0.005, () => tr.length, (v) => { tr.length = v; setDeployMax(v); }),
);
trackPane.addEventListener('toggle', render);

const hingePane = section('Flaperon hinge (pivot)',
  slider('pivot x', 0.6, 1.05, 0.005, () => pivot[0], (v) => (pivot[0] = v)),
  slider('pivot y', -0.15, 0.1, 0.005, () => pivot[1], (v) => (pivot[1] = v)),
);
hingePane.addEventListener('toggle', render);

document.getElementById('controls')!.append(
  // ── Level 1: main airfoil (foundational; fixed for now) ──
  section('① Main airfoil',
    info('airfoil', 'LS(1)-0417 (fixed)'),
    slider('blunt TE', 0.0, 0.02, 0.001, () => m.bluntThickness, (v) => (m.bluntThickness = v)),
  ),
  // ── Level 2: high-lift element design (sub-boxes are individually collapsible) ──
  section('② High-lift design',
    section('Main cove cutout',
      slider('upper cut x', 0.5, 1.0, 0.005, () => c.upperCutX, (v) => (c.upperCutX = v)),
      slider('lower cut x', 0.5, 1.0, 0.005, () => c.lowerCutX, (v) => (c.lowerCutX = v)),
      slider('cove x', 0.5, 1.0, 0.005, () => c.coveX, (v) => (c.coveX = v)),
    ),
    vanePane,
    flaperonPane,
    trackPane,
    hingePane,
  ),
  // ── Level 3: operation — deploy & flaperon angle (operating point) ──
  section('③ Operation',
    deployRow,
    flapRow,
  ),
);

// --- RANS run control (drives the CD–CL polar + flow-field overlay) ---
const runBtn = document.createElement('button');
runBtn.id = 'run-rans';
runBtn.className = 'run-rans';
runBtn.textContent = `Run RANS (α sweep: ${ALPHA_SWEEP[0]}…${ALPHA_SWEEP[ALPHA_SWEEP.length - 1]}°)`;
runBtn.onclick = runRans;
// Mesh quality for the next sweep: accurate (~48k cells) or fast preview (~7.6k).
const qualitySel = document.createElement('select');
qualitySel.className = 'rans-quality';
qualitySel.innerHTML = '<option value="accurate">accurate mesh</option><option value="fast">fast mesh</option>';
qualitySel.value = meshQuality;
qualitySel.onchange = () => { meshQuality = qualitySel.value as 'fast' | 'accurate'; };
const statusEl = document.createElement('div');
statusEl.id = 'rans-status';
statusEl.className = 'rans-status';
document.getElementById('polar-controls')!.append(runBtn, qualitySel, statusEl);

render();
drawPolar();
// Re-render whenever the canvas box actually changes size (window resize, or a layout
// reflow such as a control pane opening and bringing up the page scrollbar). The LIC is
// a separate canvas, so it must re-sync its buffer + camera or it desyncs from #view.
new ResizeObserver(() => {
  render();
  renderFlow();        // re-syncs the LIC buffer/camera (and redraws the current flow, if any)
}).observe(canvas);
