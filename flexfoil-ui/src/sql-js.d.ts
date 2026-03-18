declare module 'sql.js' {
  export interface Statement {
    run(params?: unknown[]): void;
    get(params?: unknown[]): (string | number | null | Uint8Array)[];
    getAsObject(params?: unknown[]): Record<string, string | number | null | Uint8Array>;
    step(): boolean;
    free(): void;
    bind(params?: unknown[]): boolean;
  }

  export interface Database {
    run(sql: string, params?: unknown[]): void;
    exec(sql: string, params?: unknown[]): { columns: string[]; values: (string | number | null | Uint8Array)[][] }[];
    prepare(sql: string): Statement;
    export(): Uint8Array;
    close(): void;
  }

  export interface SqlJsStatic {
    Database: new (data?: Uint8Array) => Database;
  }

  export interface InitSqlJsOptions {
    locateFile?: (file: string) => string;
  }

  export default function initSqlJs(options?: InitSqlJsOptions): Promise<SqlJsStatic>;
}
