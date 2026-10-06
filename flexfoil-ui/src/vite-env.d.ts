/// <reference types="vite/client" />

declare const __APP_VERSION__: string;

interface ImportMetaEnv {
  readonly VITE_SOLVE_RUN_ANALYTICS?: string;
  readonly VITE_AG_GRID_LICENSE_KEY: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
