import { describe, expect, it } from 'vitest';
import { parseAirfoilDat, prepareImportedAirfoil } from './airfoilImport';

const CLARK_Y = `CLARK Y AIRFOIL
 1.000000 0.000599
 0.990000 0.002969
 0.980000 0.005333
 0.970000 0.007687
 0.960000 0.010023
 0.940000 0.014624
 0.920000 0.019116
 0.900000 0.023502
 0.880000 0.027789
 0.860000 0.031974
 0.840000 0.036054
 0.820000 0.040024
 0.800000 0.043884
 0.780000 0.047628
 0.760000 0.051257
 0.740000 0.054767
 0.720000 0.058160
 0.700000 0.061433
 0.680000 0.064584
 0.660000 0.067605
 0.640000 0.070482
 0.620000 0.073206
 0.600000 0.075763
 0.580000 0.078145
 0.560000 0.080348
 0.540000 0.082371
 0.520000 0.084214
 0.500000 0.085877
 0.480000 0.087357
 0.460000 0.088643
 0.440000 0.089718
 0.420000 0.090566
 0.400000 0.091171
 0.380000 0.091521
 0.360000 0.091627
 0.340000 0.091508
 0.320000 0.091186
 0.300000 0.090680
 0.280000 0.090002
 0.260000 0.089084
 0.240000 0.087831
 0.220000 0.086143
 0.200000 0.083920
 0.180000 0.081069
 0.160000 0.077571
 0.140000 0.073436
 0.120000 0.068620
 0.100000 0.062998
 0.080000 0.056431
 0.060000 0.048757
 0.050000 0.044275
 0.040000 0.039128
 0.030000 0.033022
 0.020000 0.025374
 0.012000 0.017858
 0.008000 0.013735
 0.004000 0.008924
 0.002000 0.005803
 0.001000 0.003727
 0.000500 0.002339
 0.000000 0.000000
 0.000500 -0.004670
 0.001000 -0.005942
 0.002000 -0.007811
 0.004000 -0.010513
 0.008000 -0.014286
 0.012000 -0.016973
 0.020000 -0.020272
 0.030000 -0.022606
 0.040000 -0.024521
 0.050000 -0.026045
 0.060000 -0.027128
 0.080000 -0.028459
 0.100000 -0.029379
 0.120000 -0.029963
 0.140000 -0.030240
 0.160000 -0.030255
 0.180000 -0.030049
 0.200000 -0.029666
 0.220000 -0.029145
 0.240000 -0.028518
 0.260000 -0.027816
 0.280000 -0.027070
 0.300000 -0.026308
 0.320000 -0.025556
 0.340000 -0.024818
 0.360000 -0.024087
 0.380000 -0.023361
 0.400000 -0.022634
 0.420000 -0.021904
 0.440000 -0.021171
 0.460000 -0.020435
 0.480000 -0.019699
 0.500000 -0.018962
 0.520000 -0.018226
 0.540000 -0.017491
 0.560000 -0.016757
 0.580000 -0.016023
 0.600000 -0.015289
 0.620000 -0.014555
 0.640000 -0.013821
 0.660000 -0.013086
 0.680000 -0.012351
 0.700000 -0.011617
 0.720000 -0.010882
 0.740000 -0.010148
 0.760000 -0.009413
 0.780000 -0.008679
 0.800000 -0.007944
 0.820000 -0.007210
 0.840000 -0.006475
 0.860000 -0.005741
 0.880000 -0.005006
 0.900000 -0.004272
 0.920000 -0.003537
 0.940000 -0.002803
 0.960000 -0.002068
 0.970000 -0.001701
 0.980000 -0.001334
 0.990000 -0.000967
 1.000000 -0.000599
`;

const E205 = `E205 (10.48%)
 1.000000 0.000000
 0.996550 0.000390
 0.986490 0.001740
 0.970490 0.004270
 0.949160 0.007780
 0.922850 0.011960
 0.891750 0.016680
 0.856240 0.021990
 0.816840 0.027860
 0.774120 0.034190
 0.728660 0.040880
 0.681080 0.047770
 0.632040 0.054700
 0.582180 0.061470
 0.532170 0.067820
 0.482650 0.073420
 0.434100 0.077850
 0.386800 0.080810
 0.341010 0.082140
 0.296990 0.081770
 0.254960 0.079700
 0.215080 0.076060
 0.177640 0.071110
 0.143020 0.065070
 0.111570 0.058110
 0.083600 0.050400
 0.059370 0.042110
 0.039090 0.033440
 0.022920 0.024610
 0.010970 0.015890
 0.003310 0.007660
 0.000020 0.000550
 0.002330 -0.005060
 0.010650 -0.009880
 0.024190 -0.014200
 0.042910 -0.017760
 0.066690 -0.020530
 0.095340 -0.022520
 0.128640 -0.023780
 0.166270 -0.024360
 0.207830 -0.024350
 0.252900 -0.023840
 0.300970 -0.022920
 0.351490 -0.021680
 0.403880 -0.020210
 0.457510 -0.018590
 0.511740 -0.016890
 0.565910 -0.015160
 0.619380 -0.013450
 0.671490 -0.011800
 0.721600 -0.010230
 0.769110 -0.008760
 0.813430 -0.007400
 0.854000 -0.006140
 0.890340 -0.004970
 0.921950 -0.003800
 0.948600 -0.002520
 0.970170 -0.001250
 0.986350 -0.000360
 0.996510 -0.000030
 1.000000 0.000000
`;

const NACA_2412 = `NACA 2412
 1.000000 0.001300
 0.950000 0.011400
 0.900000 0.020800
 0.800000 0.037500
 0.700000 0.051800
 0.600000 0.063600
 0.500000 0.072400
 0.400000 0.078000
 0.300000 0.078800
 0.250000 0.076700
 0.200000 0.072600
 0.150000 0.066100
 0.100000 0.056300
 0.075000 0.049600
 0.050000 0.041300
 0.025000 0.029900
 0.012500 0.021500
 0.000000 0.000000
 0.012500 -0.016500
 0.025000 -0.022700
 0.050000 -0.030100
 0.075000 -0.034600
 0.100000 -0.037500
 0.150000 -0.041000
 0.200000 -0.042300
 0.250000 -0.042200
 0.300000 -0.041200
 0.400000 -0.038000
 0.500000 -0.033400
 0.600000 -0.027600
 0.700000 -0.021400
 0.800000 -0.015000
 0.900000 -0.008200
 0.950000 -0.004800
 1.000000 -0.001300
`;

/** Chordwise stations used to synthesise element contours, TE -> LE. */
const STATIONS = [1, 0.8, 0.6, 0.4, 0.2, 0.05, 0];

function fmt(x: number, y: number): string {
  return ` ${x.toFixed(6)} ${y.toFixed(6)}`;
}

/**
 * A closed element contour placed at `xLe` with chord `chord`.
 *
 * `startAt: 'te'` writes the usual Selig order (TE -> upper -> LE -> lower -> TE).
 * `startAt: 'nearLe'` starts one station aft of the leading edge on the lower
 * surface and wraps the whole loop, so the block starts at x < 0.1 — the shape
 * that used to be misread as a Lednicer surface.
 */
function closedElement(
  xLe: number,
  chord: number,
  { thickness = 0.08, startAt = 'te' as 'te' | 'nearLe' } = {},
): string {
  const place = (t: number) => xLe + chord * t;
  const half = (t: number) => chord * thickness * Math.sin(Math.PI * t);
  const leToTe = [...STATIONS].reverse();

  const lines =
    startAt === 'te'
      ? [
          ...STATIONS.map((t) => fmt(place(t), half(t))), // TE -> LE, upper
          ...leToTe.slice(1).map((t) => fmt(place(t), -half(t))), // LE -> TE, lower
        ]
      : [
          fmt(place(leToTe[1]), -half(leToTe[1])), // just aft of the LE, lower
          fmt(place(leToTe[0]), half(leToTe[0])), // LE
          ...leToTe.slice(1).map((t) => fmt(place(t), half(t))), // LE -> TE, upper
          ...STATIONS.slice(1, -1).map((t) => fmt(place(t), -half(t))), // TE -> back to start, lower
        ];
  return lines.join('\n');
}

/** An open Lednicer-style surface running LE -> TE. */
function lednicerSurface(sign: number): string {
  return [...STATIONS]
    .reverse()
    .map((t) => fmt(t, sign * 0.08 * Math.sin(Math.PI * t)))
    .join('\n');
}

describe('parseAirfoilDat', () => {
  it.each([
    ['Clark Y', CLARK_Y, 'clarky.dat', 'CLARK Y AIRFOIL'],
    ['E205', E205, 'e205.dat', 'E205 (10.48%)'],
    ['NACA 2412', NACA_2412, 'naca2412.dat', 'NACA 2412'],
  ])('parses real Airfoil Tools %s data', (_label, text, fileName, expectedName) => {
    const parsed = parseAirfoilDat(text, fileName);

    expect(parsed.name).toBe(expectedName);
    expect(parsed.coordinates.length).toBeGreaterThan(20);

    const leIndex = parsed.coordinates.findIndex((point) => point.x === Math.min(...parsed.coordinates.map((p) => p.x)));
    expect(leIndex).toBeGreaterThan(0);
    expect(leIndex).toBeLessThan(parsed.coordinates.length - 1);

    expect(parsed.coordinates[0].surface).toBe('upper');
    expect(parsed.coordinates[leIndex].surface).toBe('upper');
    expect(parsed.coordinates.at(-1)?.surface).toBe('lower');
    expect(parsed.coordinates[0].y).toBeGreaterThanOrEqual(-0.001);
    expect(parsed.coordinates.at(-1)?.y ?? 0).toBeLessThanOrEqual(0.001);
  });

  it('falls back to the file name when the header is missing', () => {
    const parsed = parseAirfoilDat('1.0 0.0\n0.0 0.1\n1.0 -0.0\n', 'my-foil.dat');
    expect(parsed.name).toBe('my foil');
  });

  it('handles Lednicer format (two groups separated by blank line)', () => {
    const lednicer = `CLARK X AIRFOIL
       17.       17.

 0.0000000 0.0000000
 0.0125000 0.0186000
 0.0250000 0.0276000
 0.0500000 0.0416000
 0.1000000 0.0608000
 0.2000000 0.0808000
 0.4000000 0.0900000
 0.6000000 0.0755000
 0.8000000 0.0442000
 1.0000000 0.0012000

 0.0000000 0.0000000
 0.0125000 -.0163000
 0.0250000 -.0220000
 0.0500000 -.0272000
 0.1000000 -.0310000
 0.2000000 -.0317000
 0.4000000 -.0240000
 0.6000000 -.0160000
 0.8000000 -.0080000
 1.0000000 0.0000000
`;
    const parsed = parseAirfoilDat(lednicer, 'clarkx.dat');

    expect(parsed.name).toBe('CLARK X AIRFOIL');
    // Should be converted to Selig order: TE → upper → LE → lower → TE
    expect(parsed.coordinates[0].x).toBeCloseTo(1.0);
    expect(parsed.coordinates[0].surface).toBe('upper');
    expect(parsed.coordinates.at(-1)!.x).toBeCloseTo(1.0);
    expect(parsed.coordinates.at(-1)!.surface).toBe('lower');

    const leIndex = parsed.coordinates.findIndex(
      (p) => p.x === Math.min(...parsed.coordinates.map((q) => q.x)),
    );
    expect(leIndex).toBeGreaterThan(0);
    expect(leIndex).toBeLessThan(parsed.coordinates.length - 1);
  });

  it('handles Lednicer format with decimal count line', () => {
    const lednicer = `TEST FOIL
      5.0      5.0

 0.0000 0.0000
 0.2500 0.0500
 0.5000 0.0700
 0.7500 0.0400
 1.0000 0.0000

 0.0000 0.0000
 0.2500 -0.0300
 0.5000 -0.0500
 0.7500 -0.0300
 1.0000 0.0000
`;
    const parsed = parseAirfoilDat(lednicer, 'test.dat');
    expect(parsed.name).toBe('TEST FOIL');
    // 5 upper reversed + 4 lower (skip dup LE) = 9 points
    expect(parsed.coordinates.length).toBe(9);
    expect(parsed.coordinates[0].x).toBeCloseTo(1.0);
    expect(parsed.coordinates.at(-1)!.x).toBeCloseTo(1.0);
  });
});

describe('parseAirfoilDat multi-element', () => {
  const SLAT = closedElement(-0.1, 0.15);
  const MAIN = closedElement(0, 1);
  const FLAP = closedElement(0.85, 0.45);
  const POINTS_PER_ELEMENT = STATIONS.length * 2 - 1;

  it('reports an ordinary single-element Selig file as one element', () => {
    const parsed = parseAirfoilDat(CLARK_Y, 'clarky.dat');

    expect(parsed.elements).toHaveLength(1);
    expect(parsed.coordinates).toBe(parsed.elements[0]);
  });

  it('keeps a blank-line separated slat/main/flap file as three elements', () => {
    const text = `MDA 30P-30N (mini)\n${SLAT}\n\n${MAIN}\n\n${FLAP}\n`;
    const parsed = parseAirfoilDat(text, 'threeElement.dat');

    expect(parsed.name).toBe('MDA 30P-30N (mini)');
    expect(parsed.elements).toHaveLength(3);
    expect(parsed.elements.map((e) => e.length)).toEqual([
      POINTS_PER_ELEMENT,
      POINTS_PER_ELEMENT,
      POINTS_PER_ELEMENT,
    ]);

    // Element 0 is the single-element view, not a 39-point three-loop contour.
    expect(parsed.coordinates).toBe(parsed.elements[0]);
    expect(parsed.coordinates).toHaveLength(POINTS_PER_ELEMENT);

    // Each element keeps its own chordwise station.
    const ranges = parsed.elements.map((e) => [
      Math.min(...e.map((p) => p.x)),
      Math.max(...e.map((p) => p.x)),
    ]);
    expect(ranges[0][0]).toBeCloseTo(-0.1);
    expect(ranges[0][1]).toBeCloseTo(0.05);
    expect(ranges[1][0]).toBeCloseTo(0);
    expect(ranges[1][1]).toBeCloseTo(1);
    expect(ranges[2][0]).toBeCloseTo(0.85);
    expect(ranges[2][1]).toBeCloseTo(1.3);

    // Surfaces are annotated per element, not across the whole file.
    for (const element of parsed.elements) {
      expect(element[0].surface).toBe('upper');
      expect(element.at(-1)?.surface).toBe('lower');
    }
  });

  it('keeps a comment-labelled slat/main/flap file as three elements', () => {
    // The structure of the real `30p-30n.dat` in public/airfoils: comment labels
    // are the only element boundary in the file, and each labelled block stands
    // up as a contour in its own right.
    const text = `# MDA 30P-30N (mini)\n# Slat\n${SLAT}\n# Main Element\n${MAIN}\n# Flap\n${FLAP}\n`;
    const parsed = parseAirfoilDat(text, '30p-30n.dat');

    expect(parsed.elements).toHaveLength(3);
    expect(parsed.elements.map((e) => e.length)).toEqual([
      POINTS_PER_ELEMENT,
      POINTS_PER_ELEMENT,
      POINTS_PER_ELEMENT,
    ]);
    expect(parsed.coordinates).toHaveLength(POINTS_PER_ELEMENT);
  });

  it('treats a 999.0 999.0 line as an element boundary, not a coordinate', () => {
    // The coordinate range check would have discarded the sentinel as a
    // non-coordinate in any case; what matters is that it is recognised as a
    // boundary the file declared.
    const text = `SEPARATED\n${SLAT}\n 999.0 999.0\n${MAIN}\n 999.0 999.0\n${FLAP}\n`;
    const parsed = parseAirfoilDat(text, 'separated.dat');

    expect(parsed.elements).toHaveLength(3);
    const all = parsed.elements.flat();
    expect(all).toHaveLength(3 * POINTS_PER_ELEMENT);
    expect(all.some((p) => p.x >= 999 || p.y >= 999)).toBe(false);
    expect(Math.max(...all.map((p) => p.x))).toBeCloseTo(1.3);
  });

  it('tolerates separator whitespace and format variation', () => {
    const text = `SEPARATED\n${MAIN}\n   999.   999.  \n${FLAP}\n999.000000   999.000000\n${SLAT}\n`;
    const parsed = parseAirfoilDat(text, 'separated.dat');

    expect(parsed.elements).toHaveLength(3);
    expect(parsed.elements.flat().some((p) => p.x >= 999)).toBe(false);
  });

  it('does not read a genuine two-element file as Lednicer', () => {
    // Both elements are closed loops that start near the leading edge and cover
    // the same x-range, so "both blocks start at the LE" cannot tell them from
    // Lednicer surfaces; the closure and monotonicity tests carry the decision.
    const first = closedElement(0, 1, { startAt: 'nearLe' });
    const second = closedElement(0, 1, { thickness: 0.05, startAt: 'nearLe' });
    const parsed = parseAirfoilDat(`TWO ELEMENTS\n${first}\n\n${second}\n`, 'twoElement.dat');

    // Lednicer conversion would fuse the blocks into a single contour of
    // 2 * POINTS_PER_ELEMENT - 1 points; here each block stays its own element.
    expect(parsed.elements).toHaveLength(2);
    expect(parsed.elements.map((e) => e.length)).toEqual([
      POINTS_PER_ELEMENT,
      POINTS_PER_ELEMENT,
    ]);
    expect(parsed.coordinates).toHaveLength(POINTS_PER_ELEMENT);
    // Element 0 is kept in file order rather than reversed: it still starts just
    // aft of the leading edge on the lower surface and reaches the LE next.
    expect(parsed.coordinates[0].x).toBeCloseTo(0.05);
    expect(parsed.coordinates[0].y).toBeLessThan(0);
    expect(parsed.coordinates[1].x).toBeCloseTo(0);
    // Element 1 keeps its own (thinner) section instead of being appended to 0.
    const thickness = (points: { y: number }[]) => Math.max(...points.map((p) => Math.abs(p.y)));
    expect(thickness(parsed.elements[1])).toBeLessThan(thickness(parsed.elements[0]));
  });

  it('treats comments that annotate one contour as decoration, not boundaries', () => {
    // A hand-annotated single-element file: the comments label the two surfaces
    // of one contour, so honouring them as element boundaries would report the
    // file as two elements. Each surface is longer than the minimum block length,
    // so the decision rests on the chordwise traversal test.
    const n = 20;
    const stations = Array.from({ length: n + 1 }, (_, i) => 1 - i / n); // TE -> LE
    const half = (t: number) => 0.08 * Math.sin(Math.PI * t);
    const upper = stations.map((t) => fmt(t, half(t)));
    const lower = [...stations].reverse().slice(1).map((t) => fmt(t, -half(t)));
    const text = `ANNOTATED\n# upper surface\n${upper.join('\n')}\n# lower surface\n${lower.join('\n')}\n`;

    const parsed = parseAirfoilDat(text, 'annotated.dat');

    expect(parsed.elements).toHaveLength(1);
    expect(parsed.coordinates).toHaveLength(2 * n + 1);
    expect(parsed.coordinates[0].x).toBeCloseTo(1);
    expect(parsed.coordinates.at(-1)?.x).toBeCloseTo(1);
  });

  it('still detects Lednicer format when both blocks are open LE->TE surfaces', () => {
    const text = `LEDNICER FOIL\n      7.      7.\n\n${lednicerSurface(1)}\n\n${lednicerSurface(-1)}\n`;
    const parsed = parseAirfoilDat(text, 'lednicer.dat');

    expect(parsed.name).toBe('LEDNICER FOIL');
    // Two surfaces of ONE element -> one element, in Selig order.
    expect(parsed.elements).toHaveLength(1);
    expect(parsed.coordinates).toHaveLength(STATIONS.length * 2 - 1);
    expect(parsed.coordinates[0].x).toBeCloseTo(1);
    expect(parsed.coordinates.at(-1)?.x).toBeCloseTo(1);
    expect(parsed.coordinates[0].surface).toBe('upper');
    expect(parsed.coordinates.at(-1)?.surface).toBe('lower');
  });

  it('rejoins a stray blank line inside a single element (cap21c pattern)', () => {
    // cap21c.dat in the bundled library separates its closing TE point with two
    // blank lines; the point belongs to the element and must not be dropped or
    // promoted to a one-point "element".
    const text = `CAP 21 (mini)\n${MAIN}\n\n\n 0.998900 -0.001006\n`;
    const parsed = parseAirfoilDat(text, 'cap21c.dat');

    expect(parsed.elements).toHaveLength(1);
    expect(parsed.coordinates).toHaveLength(POINTS_PER_ELEMENT + 1);
    expect(parsed.coordinates.at(-1)?.x).toBeCloseTo(0.9989);
  });

  it('reports a degenerate element declared by an explicit separator', () => {
    const text = `BROKEN\n${MAIN}\n 999.0 999.0\n 0.5 0.0\n`;

    expect(() => parseAirfoilDat(text, 'broken.dat')).toThrow(/Element 2/);
  });
});

describe('prepareImportedAirfoil', () => {
  it('uses repaneled output when a repaneler is provided', () => {
    const parsed = parseAirfoilDat(NACA_2412, 'naca2412.dat');
    const imported = prepareImportedAirfoil(parsed, 160, () => [
      { x: 1, y: 0.01 },
      { x: 0, y: 0 },
      { x: 1, y: -0.01 },
    ]);

    expect(imported.coordinates).toEqual(parsed.coordinates);
    expect(imported.panels).toEqual([
      { x: 1, y: 0.01, surface: 'upper' },
      { x: 0, y: 0, surface: 'upper' },
      { x: 1, y: -0.01, surface: 'lower' },
    ]);
  });

  it('falls back to raw coordinates when repaneling is unavailable', () => {
    const parsed = parseAirfoilDat(E205, 'e205.dat');
    const imported = prepareImportedAirfoil(parsed, 160);

    expect(imported.panels).toEqual(parsed.coordinates);
  });
});
