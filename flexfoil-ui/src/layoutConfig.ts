import type { IJsonModel } from 'flexlayout-react';

export const PANELS = [
  { id: 'canvas', name: 'Airfoil Canvas', component: 'canvas' },
  { id: 'library', name: 'Airfoil Library', component: 'library' },
  { id: 'control', name: 'Geometry Control', component: 'control' },
  { id: 'spacing', name: 'Spacing', component: 'spacing' },
  { id: 'properties', name: 'Properties', component: 'properties' },
  { id: 'solve', name: 'Solve', component: 'solve' },
  { id: 'polar', name: 'Polar Plot', component: 'polar' },
  { id: 'visualization', name: 'Visualization', component: 'visualization' },
  { id: 'data-explorer', name: 'Data Explorer', component: 'data-explorer' },
  { id: 'plot-builder', name: 'Plot Builder', component: 'plot-builder' },
  { id: 'case-logs', name: 'Case Logs', component: 'case-logs' },
] as const;

export type PanelId = (typeof PANELS)[number]['id'];

export const defaultLayoutJson: IJsonModel = {
  global: {
    tabSetMinWidth: 200,
    tabSetMinHeight: 150,
    borderMinSize: 100,
    splitterSize: 6,
  },
  borders: [],
  layout: {
    type: 'row',
    weight: 100,
    children: [
      {
        type: 'row',
        weight: 65,
        children: [
          {
            type: 'tabset',
            weight: 70,
            children: [
              { type: 'tab', id: 'canvas', name: 'Airfoil Canvas', component: 'canvas' },
            ],
          },
          {
            type: 'tabset',
            weight: 30,
            children: [
              { type: 'tab', id: 'polar', name: 'Polar Plot', component: 'polar' },
              { type: 'tab', id: 'plot-builder', name: 'Plot Builder', component: 'plot-builder' },
              { type: 'tab', id: 'data-explorer', name: 'Data Explorer', component: 'data-explorer' },
              { type: 'tab', id: 'case-logs', name: 'Case Logs', component: 'case-logs' },
            ],
          },
        ],
      },
      {
        type: 'row',
        weight: 35,
        children: [
          {
            type: 'tabset',
            weight: 50,
            children: [
              { type: 'tab', id: 'solve', name: 'Solve', component: 'solve' },
              { type: 'tab', id: 'properties', name: 'Properties', component: 'properties' },
            ],
          },
          {
            type: 'tabset',
            weight: 50,
            children: [
              { type: 'tab', id: 'library', name: 'Airfoil Library', component: 'library' },
              { type: 'tab', id: 'control', name: 'Geometry Control', component: 'control' },
              { type: 'tab', id: 'visualization', name: 'Visualization', component: 'visualization' },
              { type: 'tab', id: 'spacing', name: 'Spacing', component: 'spacing' },
            ],
          },
        ],
      },
    ],
  },
};

export function isPanelId(value: string | null | undefined): value is PanelId {
  return !!value && PANELS.some((panel) => panel.id === value);
}
