import type { RunRow } from '../types';

export type PlotFieldKind = 'numeric' | 'categorical' | 'boolean' | 'temporal';

export interface PlotFieldMeta {
  key: keyof RunRow;
  label: string;
  kind: PlotFieldKind;
  autoGroupInvariant?: boolean;
}

export const PLOT_FIELDS: PlotFieldMeta[] = [
  { key: 'id', label: 'ID', kind: 'numeric' },
  { key: 'airfoil_name', label: 'Airfoil', kind: 'categorical', autoGroupInvariant: true },
  { key: 'airfoil_hash', label: 'Airfoil Hash', kind: 'categorical', autoGroupInvariant: true },
  { key: 'alpha', label: 'Alpha', kind: 'numeric' },
  { key: 'reynolds', label: 'Re', kind: 'numeric', autoGroupInvariant: true },
  { key: 'mach', label: 'Mach', kind: 'numeric', autoGroupInvariant: true },
  { key: 'ncrit', label: 'Ncrit', kind: 'numeric', autoGroupInvariant: true },
  { key: 'n_panels', label: 'N Panels', kind: 'numeric', autoGroupInvariant: true },
  { key: 'max_iter', label: 'Max Iter', kind: 'numeric', autoGroupInvariant: true },
  { key: 'cl', label: 'CL', kind: 'numeric' },
  { key: 'cd', label: 'CD', kind: 'numeric' },
  { key: 'cm', label: 'CM', kind: 'numeric' },
  { key: 'ld', label: 'L/D', kind: 'numeric' },
  { key: 'converged', label: 'Converged', kind: 'boolean' },
  { key: 'iterations', label: 'Iterations', kind: 'numeric' },
  { key: 'residual', label: 'Residual', kind: 'numeric' },
  { key: 'x_tr_upper', label: 'Xtr Upper', kind: 'numeric' },
  { key: 'x_tr_lower', label: 'Xtr Lower', kind: 'numeric' },
  { key: 'flap_deflection', label: 'Flap δ', kind: 'numeric' },
  { key: 'flap_hinge_x', label: 'Flap x/c', kind: 'numeric', autoGroupInvariant: true },
  { key: 'solver_mode', label: 'Solver', kind: 'categorical', autoGroupInvariant: true },
  { key: 'created_at', label: 'Created', kind: 'temporal' },
  { key: 'session_id', label: 'Session', kind: 'categorical' },
];

export const NUMERIC_PLOT_FIELDS = PLOT_FIELDS.filter(
  (field) => field.kind === 'numeric',
);

const DISCRETE_PLOT_FIELDS = PLOT_FIELDS.filter(
  (field) => field.kind !== 'numeric',
);

export const ENCODING_PLOT_FIELDS = PLOT_FIELDS;

export const AUTO_GROUP_INVARIANT_FIELDS = PLOT_FIELDS.filter(
  (field) => field.autoGroupInvariant,
).map((field) => field.key);

export function getPlotFieldMeta(key: keyof RunRow): PlotFieldMeta | undefined {
  return PLOT_FIELDS.find((field) => field.key === key);
}

export function getPlotFieldLabel(key: keyof RunRow): string {
  return getPlotFieldMeta(key)?.label ?? String(key);
}

export function isNumericPlotField(key: keyof RunRow): boolean {
  return getPlotFieldMeta(key)?.kind === 'numeric';
}

function isDiscretePlotField(key: keyof RunRow): boolean {
  return !isNumericPlotField(key);
}
