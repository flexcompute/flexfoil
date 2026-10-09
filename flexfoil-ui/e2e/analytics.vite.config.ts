import {defineConfig} from 'vite';
import react from '@vitejs/plugin-react';
import {fileURLToPath} from 'node:url';
const root = fileURLToPath(new URL('../../', import.meta.url));
const wasmFixture = fileURLToPath(new URL('../../e2e/analytics/wasm-fixture.ts', import.meta.url));
export default defineConfig({root, plugins: [react(), {
  name: 'offline-solver-for-analytics', enforce: 'pre',
  resolveId(source) {if (source.endsWith('/lib/wasm') || source === './wasm') return wasmFixture;},
}], resolve: {dedupe: ['react', 'react-dom'], alias: {
  react: fileURLToPath(new URL('../node_modules/react', import.meta.url)),
  'react-dom': fileURLToPath(new URL('../node_modules/react-dom', import.meta.url)),
}}});
