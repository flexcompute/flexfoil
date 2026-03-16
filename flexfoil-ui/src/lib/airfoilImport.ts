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

function isLikelyCountLine(line: string): boolean {
  const trimmed = line.trim();
  if (!/^\d+\s+\d+$/.test(trimmed)) return false;

  const [a, b] = trimmed.split(/\s+/).map(Number);
  return a > 2 && b > 2;
}

function parseCoordinateLine(line: string): AirfoilPoint | null {
  const trimmed = line.trim();
  if (!trimmed || isLikelyCountLine(trimmed)) {
    return null;
  }

  const parts = trimmed.replace(/,/g, ' ').split(/\s+/);
  if (parts.length < 2) {
    return null;
  }

  const x = Number(parts[0]);
  const y = Number(parts[1]);
  if (!Number.isFinite(x) || !Number.isFinite(y)) {
    return null;
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

export function parseAirfoilDat(text: string, fileName: string): ParsedAirfoilFile {
  const lines = text.split(/\r?\n/);
  const coordinates: AirfoilPoint[] = [];
  let header: string | null = null;

  for (const line of lines) {
    const coordinate = parseCoordinateLine(line);
    if (coordinate) {
      coordinates.push(coordinate);
      continue;
    }

    const trimmed = line.trim();
    if (!trimmed || header !== null || coordinates.length > 0) {
      continue;
    }

    header = trimmed;
  }

  if (coordinates.length < 3) {
    throw new Error('Expected at least 3 airfoil coordinates.');
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
    if (repaneled.length > 0) {
      panels = annotateSurfaces(repaneled);
    }
  }

  return {
    ...parsed,
    panels,
  };
}
