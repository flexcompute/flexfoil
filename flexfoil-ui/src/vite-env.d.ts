/// <reference types="vite/client" />

interface ImportMetaEnv {
  readonly VITE_AG_GRID_LICENSE_KEY: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
