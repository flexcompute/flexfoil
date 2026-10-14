// Exercise the shipping WASM exports against the same live reference and native results.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const wasm = require(path.resolve(process.argv[2]));
const report = JSON.parse(fs.readFileSync(process.argv[3], 'utf8'));
assert.equal(report.cases.length, Number(process.argv[4]));
const tolerances = { cl: 2e-4, cd: 1e-5, cm: 1e-4, x_tr_upper: 2e-4, x_tr_lower: 2e-4 };
function close(actual, expected, tolerance, label) {
  assert.ok(Number.isFinite(actual) && Number.isFinite(expected) && Math.abs(actual - expected) <= tolerance,
    `${label}: ${actual} versus ${expected}, tolerance ${tolerance}`);
}
const failures = [];
for (const c of report.cases) {
  try {
    const result = wasm.analyze_airfoil_faithful(new Float64Array(c.coords), c.alpha, c.reynolds,
      c.mach, c.ncrit, 100, 1, ...c.trips);
    assert.equal(result.success, true);
    assert.equal(result.converged, true);
    for (const [key, tolerance] of Object.entries(tolerances)) {
      close(result[key], c.reference[key], tolerance, `${c.name} WASM/reference ${key}`);
      close(result[key], c.rust[key], 1e-10, `${c.name} WASM/native ${key}`);
    }
    assert.equal(result.cp.length, c.cp_reference.length);
    assert.equal(result.cp_x.length, result.cp.length);
    for (let i = 0; i < result.cp.length; i++) {
      close(result.cp_x[i], c.cp_reference[i][0], 5.1e-6, `${c.name} Cp location ${i}`);
      close(result.cp[i], c.cp_reference[i][1], 2e-4, `${c.name} Cp ${i}`);
    }
  } catch (error) {
    failures.push(`${c.name}: ${error.message}`);
  }
}
assert.deepEqual(failures, []);
console.log(`WASM/reference/native: ${report.cases.length} complete cases passed`);
