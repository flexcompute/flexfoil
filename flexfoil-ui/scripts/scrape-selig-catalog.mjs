#!/usr/bin/env node
/**
 * One-time scraper: fetches the UIUC Airfoil Coordinates Database HTML page
 * and extracts every .dat link + its description into a JSON catalog file.
 *
 * Usage:  node scripts/scrape-selig-catalog.mjs
 * Output: src/data/seligCatalog.json
 */

import { writeFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const OUTPUT = join(__dirname, '..', 'src', 'data', 'seligCatalog.json');

const UIUC_URL = 'https://m-selig.ae.illinois.edu/ads/coord_database.html';

async function main() {
  console.log(`Fetching ${UIUC_URL} ...`);
  const res = await fetch(UIUC_URL);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const html = await res.text();

  // HTML pattern (relative links):
  //   <a href="coord/FILENAME.dat">FILENAME.dat</a> \ DESCRIPTION \ <a href="afplots/...
  //   <a href="coord_updates/FILENAME.dat">FILENAME.dat</a> \ DESCRIPTION ...
  const seen = new Set();
  const catalog = [];

  // Match: <a href="coord/FILE.dat">...</a> followed by backslash-separated description
  const re = /href="coord(?:_updates)?\/([^"]+\.dat)"[^>]*>[^<]+<\/a>\s*\\\s*([^\\<]*)/gi;

  let m;
  while ((m = re.exec(html)) !== null) {
    const file = m[1].replace(/\.dat$/i, '');
    const rawDesc = m[2].replace(/\s+/g, ' ').trim();

    if (seen.has(file)) continue;
    seen.add(file);

    catalog.push({ file, name: rawDesc || file });
  }

  // Sort alphabetically by file
  catalog.sort((a, b) => a.file.localeCompare(b.file));

  writeFileSync(OUTPUT, JSON.stringify(catalog, null, 2) + '\n');
  console.log(`Wrote ${catalog.length} airfoils to ${OUTPUT}`);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
