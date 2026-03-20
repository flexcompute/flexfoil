import { parseAirfoilDat, type ParsedAirfoilFile } from './airfoilImport';

export interface SeligCatalogEntry {
  file: string;
  name: string;
}

let _catalog: SeligCatalogEntry[] | null = null;
let _catalogPromise: Promise<SeligCatalogEntry[]> | null = null;

export function getSeligCatalog(): Promise<SeligCatalogEntry[]> {
  if (_catalog) return Promise.resolve(_catalog);
  if (_catalogPromise) return _catalogPromise;

  _catalogPromise = import('../data/seligCatalog.json').then((mod) => {
    _catalog = mod.default as SeligCatalogEntry[];
    return _catalog;
  });
  return _catalogPromise;
}

const _datCache = new Map<string, string>();

async function fetchDatText(file: string): Promise<string> {
  const cached = _datCache.get(file);
  if (cached) return cached;

  const base = import.meta.env.BASE_URL ?? '/';
  const resp = await fetch(`${base}airfoils/${file}.dat`);
  if (!resp.ok) {
    throw new Error(`Failed to fetch ${file}.dat (HTTP ${resp.status})`);
  }
  const text = await resp.text();
  _datCache.set(file, text);
  return text;
}

export async function fetchSeligAirfoil(file: string): Promise<ParsedAirfoilFile> {
  const text = await fetchDatText(file);
  return parseAirfoilDat(text, `${file}.dat`);
}

export function getRandomCatalogEntry(catalog: SeligCatalogEntry[]): SeligCatalogEntry {
  return catalog[Math.floor(Math.random() * catalog.length)];
}
