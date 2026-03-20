import type { PanelId } from '../layoutConfig';
import type { ControlMode } from '../types';
import type { TourId } from '../onboarding';

export type SearchCategory = 'panel' | 'feature' | 'action' | 'tour';

export interface SearchItem {
  id: string;
  label: string;
  description: string;
  category: SearchCategory;
  keywords: string[];
  panelId?: PanelId;
  /** Extra callback after the panel is opened (e.g. set control mode). */
  postAction?: () => void;
  /** External URL to open instead of a panel. */
  href?: string;
  shortcut?: string;
}

export function buildSearchIndex(deps: {
  setControlMode: (mode: ControlMode) => void;
  startTour: (id: TourId, force?: boolean) => void;
  toggleTheme: () => void;
  undo: () => void;
  redo: () => void;
  resetLayout: () => void;
}): SearchItem[] {
  const { setControlMode, startTour, toggleTheme, undo, redo, resetLayout } = deps;

  return [
    // ── Panels ──────────────────────────────────────────
    {
      id: 'panel:canvas',
      label: 'Airfoil Canvas',
      description: 'Interactive airfoil shape viewer',
      category: 'panel',
      keywords: ['airfoil', 'shape', 'profile', 'geometry', 'wing', 'view', 'draw', 'canvas', 'foil'],
      panelId: 'canvas',
    },
    {
      id: 'panel:library',
      label: 'Airfoil Library (NACA / Selig)',
      description: 'Browse and select NACA / Selig airfoils',
      category: 'panel',
      keywords: ['NACA', 'library', 'browse', 'select', 'database', 'Selig', 'catalog', 'airfoil', 'choose', 'pick', 'NACA 4-digit'],
      panelId: 'library',
    },
    {
      id: 'panel:control',
      label: 'Geometry Control',
      description: 'Edit airfoil shape, camber, thickness, flaps',
      category: 'panel',
      keywords: ['edit', 'modify', 'geometry', 'control', 'shape', 'design'],
      panelId: 'control',
    },
    {
      id: 'panel:spacing',
      label: 'Spacing',
      description: 'Panel distribution and mesh refinement',
      category: 'panel',
      keywords: ['spacing', 'mesh', 'panel', 'SSP', 'grid', 'resolution', 'refinement', 'knots', 'distribution', 'points'],
      panelId: 'spacing',
    },
    {
      id: 'panel:properties',
      label: 'Properties',
      description: 'Reynolds, Mach, Ncrit, solver settings',
      category: 'panel',
      keywords: ['Reynolds', 'Mach', 'ncrit', 'solver', 'settings', 'viscous', 'inviscid', 'parameters', 'flight conditions', 'Re', 'properties', 'configuration'],
      panelId: 'properties',
    },
    {
      id: 'panel:solve',
      label: 'Solve',
      description: 'Run single-point or polar sweep analysis',
      category: 'panel',
      keywords: ['run', 'solve', 'compute', 'analyze', 'alpha', 'angle of attack', 'sweep', 'converge', 'calculate', 'simulation'],
      panelId: 'solve',
    },
    {
      id: 'panel:polar',
      label: 'Polar Plot',
      description: 'Lift, drag, and moment polar curves',
      category: 'panel',
      keywords: ['polar', 'polars', 'lift', 'drag', 'L/D', 'Cl', 'Cd', 'Cm', 'chart', 'performance', 'plot', 'aerodynamic', 'coefficient', 'polar diagram'],
      panelId: 'polar',
    },
    {
      id: 'panel:visualization',
      label: 'Visualization',
      description: 'Flow field, streamlines, Cp, boundary layer',
      category: 'panel',
      keywords: ['flow', 'streamlines', 'smoke', 'pressure', 'Cp', 'forces', 'boundary layer', 'wake', 'displacement thickness', 'visualization', 'field', 'contours'],
      panelId: 'visualization',
    },
    {
      id: 'panel:data-explorer',
      label: 'Data Explorer',
      description: 'Tabular data, filters, SPLOM correlogram',
      category: 'panel',
      keywords: ['table', 'data', 'filter', 'sort', 'correlogram', 'SPLOM', 'scatter matrix', 'export', 'explorer', 'grid', 'runs'],
      panelId: 'data-explorer',
    },
    {
      id: 'panel:plot-builder',
      label: 'Plot Builder',
      description: 'Custom scatter, line, bar, histogram charts',
      category: 'panel',
      keywords: ['chart', 'graph', 'scatter', 'line', 'bar', 'histogram', 'custom plot', 'plot builder', 'x-y'],
      panelId: 'plot-builder',
    },
    {
      id: 'panel:case-logs',
      label: 'Case Logs',
      description: 'Solver convergence history and output',
      category: 'panel',
      keywords: ['logs', 'history', 'solver output', 'convergence', 'debug', 'console', 'case'],
      panelId: 'case-logs',
    },
    {
      id: 'panel:distributions',
      label: 'Distributions',
      description: 'Surface distribution plots (Cp, Cf, BL quantities) vs x/c, y, or arc-length',
      category: 'panel',
      keywords: ['cp', 'cf', 'distribution', 'surface', 'boundary layer', 'pressure', 'friction', 'thickness', 'shape factor'],
      panelId: 'distributions',
    },

    // ── Sub-features (open panel + set mode) ────────────
    {
      id: 'feature:flap-design',
      label: 'Flap Design',
      description: 'Design control surfaces — hinge, deflection, TE gap',
      category: 'feature',
      keywords: ['flap', 'control surface', 'hinge', 'deflection', 'aileron', 'elevator', 'rudder', 'trailing edge', 'TE gap', 'fold', 'GDES', 'geometry design'],
      panelId: 'control',
      postAction: () => setControlMode('geometry-design'),
    },
    {
      id: 'feature:inverse-design',
      label: 'Inverse Design (QDES)',
      description: 'Specify target pressure distribution',
      category: 'feature',
      keywords: ['inverse', 'QDES', 'pressure distribution', 'target Cp', 'inverse design', 'pressure', 'specify'],
      panelId: 'control',
      postAction: () => setControlMode('inverse-design'),
    },
    {
      id: 'feature:camber-spline',
      label: 'Camber Spline Editor',
      description: 'Edit camber line via control points',
      category: 'feature',
      keywords: ['camber', 'spline', 'meanline', 'curvature', 'camber line', 'B-spline'],
      panelId: 'control',
      postAction: () => setControlMode('camber-spline'),
    },
    {
      id: 'feature:thickness-spline',
      label: 'Thickness Spline Editor',
      description: 'Edit thickness distribution via control points',
      category: 'feature',
      keywords: ['thickness', 'spline', 'distribution', 'B-spline', 'profile'],
      panelId: 'control',
      postAction: () => setControlMode('thickness-spline'),
    },
    {
      id: 'feature:parameters',
      label: 'Parameter Controls',
      description: 'NACA parameters, camber/thickness scaling',
      category: 'feature',
      keywords: ['parameters', 'NACA', 'camber scale', 'thickness scale', 'sliders'],
      panelId: 'control',
      postAction: () => setControlMode('parameters'),
    },
    {
      id: 'feature:polar-sweep',
      label: 'Polar Sweep',
      description: 'Sweep angle of attack for a drag polar',
      category: 'feature',
      keywords: ['polar sweep', 'alpha sweep', 'angle of attack sweep', 'drag polar', 'batch', 'range'],
      panelId: 'solve',
    },
    {
      id: 'feature:multi-sweep',
      label: 'Multi-Parameter Sweep',
      description: 'Sweep Re, Mach, Ncrit, flap deflection, or hinge x/c',
      category: 'feature',
      keywords: ['multi sweep', 'parameter sweep', 'Reynolds sweep', 'Mach sweep', 'matrix sweep', 'grid sweep', 'batch'],
      panelId: 'solve',
    },

    // ── Actions ─────────────────────────────────────────
    {
      id: 'action:undo',
      label: 'Undo',
      description: 'Undo the last change',
      category: 'action',
      keywords: ['undo', 'revert', 'back'],
      shortcut: 'Cmd+Z',
      postAction: undo,
    },
    {
      id: 'action:redo',
      label: 'Redo',
      description: 'Redo the last undone change',
      category: 'action',
      keywords: ['redo', 'forward'],
      shortcut: 'Cmd+Shift+Z',
      postAction: redo,
    },
    {
      id: 'action:toggle-theme',
      label: 'Toggle Dark / Light Mode',
      description: 'Switch between dark and light themes',
      category: 'action',
      keywords: ['dark mode', 'light mode', 'theme', 'toggle', 'appearance', 'night'],
      postAction: toggleTheme,
    },
    {
      id: 'action:reset-layout',
      label: 'Reset Layout',
      description: 'Restore the default panel arrangement',
      category: 'action',
      keywords: ['reset', 'layout', 'default', 'restore', 'arrangement'],
      postAction: resetLayout,
    },
    {
      id: 'action:docs',
      label: 'Open Documentation',
      description: 'FlexFoil docs site',
      category: 'action',
      keywords: ['docs', 'documentation', 'help', 'manual', 'guide', 'reference'],
      href: 'https://foil.flexcompute.com/docs/',
    },

    // ── Tours ───────────────────────────────────────────
    {
      id: 'tour:welcome',
      label: 'Welcome Tour',
      description: 'Guided introduction to FlexFoil',
      category: 'tour',
      keywords: ['welcome', 'tour', 'tutorial', 'onboarding', 'getting started', 'intro'],
      postAction: () => startTour('welcome', true),
    },
    {
      id: 'tour:airfoil-editing',
      label: 'Airfoil Editing Guide',
      description: 'Learn to edit airfoil geometry',
      category: 'tour',
      keywords: ['editing', 'guide', 'tutorial', 'airfoil editing'],
      postAction: () => startTour('airfoilEditing', true),
    },
    {
      id: 'tour:solving',
      label: 'Solving Guide',
      description: 'Learn to run analyses and polar sweeps',
      category: 'tour',
      keywords: ['solving', 'guide', 'tutorial', 'analysis guide'],
      postAction: () => startTour('solving', true),
    },
    {
      id: 'tour:data-explorer',
      label: 'Data Explorer Guide',
      description: 'Learn to explore and plot run data',
      category: 'tour',
      keywords: ['data', 'guide', 'tutorial', 'explorer guide'],
      postAction: () => startTour('dataExplorer', true),
    },
  ];
}
