import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import wasm from 'vite-plugin-wasm'
import topLevelAwait from 'vite-plugin-top-level-await'
import { readFileSync } from 'fs'

const pkg = JSON.parse(readFileSync('./package.json', 'utf-8'))

// https://vite.dev/config/
export default defineConfig({
  base: process.env.BASE_URL || '/',
  define: {
    __APP_VERSION__: JSON.stringify(pkg.version),
  },
  plugins: [
    react(),
    wasm(),
    topLevelAwait(),
  ],
  optimizeDeps: {
    exclude: ['rustfoil-wasm'],
  },
  server: {
    fs: {
      allow: ['..'],
    },
    proxy: {
      '/api/uiuc-proxy': {
        target: 'https://m-selig.ae.illinois.edu',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api\/uiuc-proxy/, '/ads/coord'),
      },
    },
  },
})
