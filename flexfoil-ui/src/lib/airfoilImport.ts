import type { AirfoilPoint } from '../types';

export interface ParsedAirfoilFile {
  name: string;
  coordinates: AirfoilPoint[];
}

export interface ImportedAirfoilState extends ParsedAirfoilFile {
  panels: AirfoilPoint[];
}

type RepanelAirfoil = (coordinates: AirfoilPoint[], nPanels: number) => AirfoilPoint[];

function fallbackAirfoilName(fileName: string): string {
  return fileName
    .replace(/\.[^.]+$/, '')
    .replace(/[_-]+/g, ' ')
    .trim() || 'Imported Airfoil';
}

/** Matches Lednicer-style count lines like "61 61", "61.0  61.0", "17. 17." */
function isLikelyCountLine(line: string): boolean {
  const trimmed = line.trim();
  if (!/^\d+\.?\d*\s+\d+\.?\d*$/.test(trimmed)) return false;

  const [a, b] = trimmed.split(/\s+/).map(Number);
  return a > 2 && b > 2;
}

function parseCoordinateLine(line: string, lineNumber: number, allowHeader: boolean): AirfoilPoint | null {
  const trimmed = line.trim();
  if (!trimmed || isLikelyCountLine(trimmed)) {
    return null;
  }

  const parts = trimmed.replace(/,/g, ' ').split(/\s+/);
  const numericLike = /^[+-]?(?:\d|\.\d|NaN\b|Inf(?:inity)?\b)/i.test(parts[0]);
  const nonFiniteToken = /^[+-]?(?:NaN|Inf(?:inity)?)$/i.test(parts[0]);
  if (parts.length < 2 || parts[1] === '') {
    if (nonFiniteToken || (!allowHeader && numericLike)) {
      throw new Error(`Missing coordinate at line ${lineNumber}.`);
    }
    return null;
  }

  const x = Number(parts[0]);
  const y = Number(parts[1]);
  if (!Number.isFinite(x) || !Number.isFinite(y)) {
    if (nonFiniteToken || (!allowHeader && numericLike)) {
      throw new Error(`Invalid or non-finite coordinate at line ${lineNumber}.`);
    }
    return null;
  }

  if (x > 1.5 || x < -0.5 || y > 1.5 || y < -1.5) {
    throw new Error(`Coordinate outside supported chord scale at line ${lineNumber}; normalize the airfoil before importing.`);
  }

  return { x, y };
}

function findLeadingEdgeIndex(coordinates: AirfoilPoint[]): number {
  let leIndex = 0;
  for (let i = 1; i < coordinates.length; i += 1) {
    if (coordinates[i].x < coordinates[leIndex].x) {
      leIndex = i;
    }
  }
  return leIndex;
}

function annotateSurfaces(coordinates: AirfoilPoint[]): AirfoilPoint[] {
  const leIndex = findLeadingEdgeIndex(coordinates);

  return coordinates.map((point, index) => ({
    ...point,
    surface: index <= leIndex ? 'upper' : 'lower',
  }));
}

/**
 * Parse coordinate groups separated by blank lines / count lines.
 * Returns { header, groups } where each group is an array of AirfoilPoints.
 */
function parseGroups(text: string): { header: string | null; groups: AirfoilPoint[][] } {
  const lines = text.split(/\r?\n/);
  let header: string | null = null;
  const groups: AirfoilPoint[][] = [];
  let current: AirfoilPoint[] = [];

  for (const [index, line] of lines.entries()) {
    if (/^(#|!|\/\/)/.test(line.trim())) continue;
    const coord = parseCoordinateLine(line, index + 1, header === null && groups.length === 0 && current.length === 0);
    if (coord) {
      current.push(coord);
      continue;
    }

    // Non-coordinate line — flush the current group if non-empty
    if (current.length > 0) {
      groups.push(current);
      current = [];
    }

    const trimmed = line.trim();
    if (!trimmed || isLikelyCountLine(trimmed)) continue;

    // First non-blank, non-count, non-coordinate line is the header
    if (header === null && groups.length === 0) {
      header = trimmed;
    }
  }
  if (current.length > 0) groups.push(current);

  return { header, groups };
}

/**
 * Detect whether coordinate groups represent Lednicer format (two groups,
 * upper LE→TE then lower LE→TE) and convert to Selig order if so.
 */
function lednicerToSelig(groups: AirfoilPoint[][]): AirfoilPoint[] | null {
  if (groups.length !== 2) return null;

  const [upper, lower] = groups;
  if (upper.length < 2 || lower.length < 2) return null;

  // Lednicer: both groups start near LE (x ≈ 0) and end near TE (x ≈ 1)
  const upperStartsAtLE = upper[0].x < 0.1;
  const lowerStartsAtLE = lower[0].x < 0.1;
  if (!upperStartsAtLE || !lowerStartsAtLE) return null;

  // Reverse upper (TE→LE) and append lower (LE→TE), skipping duplicate LE point
  const reversed = [...upper].reverse();
  return [...reversed, ...lower.slice(1)];
}

export function parseAirfoilDat(text: string, fileName: string): ParsedAirfoilFile {
  const { header, groups } = parseGroups(text);

  let coordinates: AirfoilPoint[];

  // Try Lednicer (two-group) conversion first
  const fromLednicer = lednicerToSelig(groups);
  if (fromLednicer) {
    coordinates = fromLednicer;
  } else {
    // Selig format — flatten all groups into one list
    coordinates = groups.flat();
  }

  if (coordinates.length < 3) {
    throw new Error('Expected at least 3 airfoil coordinates.');
  }

  for (let i = 1; i < coordinates.length; i += 1) {
    const previous = coordinates[i - 1];
    if (coordinates[i].x === previous.x && coordinates[i].y === previous.y) {
      throw new Error(`Duplicate adjacent coordinates at points ${i} and ${i + 1}; remove one point before importing.`);
    }
  }

  const leIndex = findLeadingEdgeIndex(coordinates);
  if (leIndex === 0 || leIndex === coordinates.length - 1) {
    throw new Error('Expected a full airfoil loop ordered TE -> upper -> LE -> lower -> TE.');
  }

  return {
    name: header ?? fallbackAirfoilName(fileName),
    coordinates: annotateSurfaces(coordinates),
  };
}

export function prepareImportedAirfoil(
  parsed: ParsedAirfoilFile,
  nPanels: number,
  repanelAirfoil?: RepanelAirfoil,
): ImportedAirfoilState {
  let panels = parsed.coordinates;

  if (repanelAirfoil) {
    const repaneled = repanelAirfoil(parsed.coordinates, nPanels);
    if (repaneled.length < 3 || repaneled.some((p) => !Number.isFinite(p.x) || !Number.isFinite(p.y))) {
      throw new Error('Repaneling failed; original airfoil has been kept. Check the input coordinates and panel spacing.');
    }
    panels = annotateSurfaces(repaneled);
  }

  return {
    ...parsed,
    panels,
  };
}
