/**
 * Standalone LIC debug harness (served at /lic-debug.html). Loads the most recent
 * RANS flow field from the bridge's replay endpoint (GET /api/rans/last-flowfield)
 * and renders it with the same LicView the design tool uses — so the LIC pipeline
 * can be debugged with a plain browser refresh, no RANS re-solve.
 *
 *   "Raw colors" → render the mesh flat-shaded by Mach (isolates geometry/camera).
 *   "LIC"        → the full two-pass render.
 *   sliders      → tune streak length (passStep) / noise frequency (texRepeat).
 */
import { LicView } from './lic';
import type { FlowMesh } from './lic';

const VIEW = { xmin: -0.1, xmax: 1.2, ymin: -0.4, ymax: 0.2 };
const W = 900;
const H = 360;

const canvas = document.getElementById('lic-canvas') as HTMLCanvasElement;
const view = new LicView(canvas);
let current: FlowMesh | null = null;

function log(msg: string): void {
  document.getElementById('log')!.textContent = msg;
  console.log('[lic-debug]', msg);
}

/** Match the design-tool letterbox so the field fills the box the same way. */
function setupCamera(): void {
  const dpr = window.devicePixelRatio || 1;
  canvas.style.width = W + 'px';
  canvas.style.height = H + 'px';
  view.resize(Math.round(W * dpr), Math.round(H * dpr));
  const s = Math.min(W / (VIEW.xmax - VIEW.xmin), H / (VIEW.ymax - VIEW.ymin));
  const ox = (W - s * (VIEW.xmax - VIEW.xmin)) / 2;
  const oy = (H - s * (VIEW.ymax - VIEW.ymin)) / 2;
  view.setCamera(VIEW.xmin - ox / s, VIEW.xmax + ox / s, VIEW.ymin - oy / s, VIEW.ymax + oy / s);
}

async function load(): Promise<void> {
  log('loading /api/rans/last-flowfield …');
  try {
    const r = await fetch('/api/rans/last-flowfield').then((res) => res.json());
    if (r.error) { log('error: ' + r.error); return; }
    current = r.flowField as FlowMesh;
    setupCamera();
    view.setFlow(current);
    view.render();
    const b = current.bounds as number[];
    log(`loaded ${r.source}\n${current.points.length / 2} pts, ${current.tris.length / 3} tris, `
      + `mach [${current.machRange.map((x) => x.toFixed(3))}]\nbounds [${b.map((x) => x.toFixed(2))}]`);
  } catch (e) {
    log('fetch failed: ' + (e as Error).message + '  (is rans_server.py running with the GET endpoint?)');
  }
}

document.getElementById('reload')!.onclick = load;
document.getElementById('colors')!.onclick = () => { if (current) { setupCamera(); view.setFlow(current); view.debugColors(); } };
document.getElementById('lic')!.onclick = () => { if (current) { setupCamera(); view.setFlow(current); view.render(); } };

const ps = document.getElementById('passStep') as HTMLInputElement;
const tr = document.getElementById('texRepeat') as HTMLInputElement;
ps.oninput = () => {
  document.getElementById('passStepV')!.textContent = ps.value;
  view.setParams({ passStep: +ps.value });
  view.render();
};
tr.oninput = () => {
  document.getElementById('texRepeatV')!.textContent = tr.value;
  view.setParams({ textureRepeat: +tr.value });
  view.render();   // render() re-runs pass 1, baking the new noise frequency
};

window.addEventListener('resize', () => { if (current) { setupCamera(); view.render(); } });
load();
