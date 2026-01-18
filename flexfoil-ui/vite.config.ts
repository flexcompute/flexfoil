import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import wasm from 'vite-plugin-wasm'
import topLevelAwait from 'vite-plugin-top-level-await'
import path from 'path'

// https://vite.dev/config/
export default defineConfig({
  base: process.env.BASE_URL || '/',
  plugins: [
    react(),
    wasm(),
    topLevelAwait(),
  ],
  resolve: {
    alias: {
      'rustfoil-wasm': path.resolve(__dirname, 'src/lib/wasm/rustfoil_wasm.js'),
    },
  },
  optimizeDeps: {
    exclude: ['rustfoil-wasm'],
  },
  server: {
    fs: {
      // Allow serving files from the pkg directory
      allow: ['..'],
    },
  },
})
