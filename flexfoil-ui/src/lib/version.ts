declare const __APP_VERSION__: string;

export const APP_VERSION: string =
  typeof __APP_VERSION__ !== 'undefined' ? __APP_VERSION__ : '0.0.0-dev';

export type ChangeCategory = 'added' | 'changed' | 'fixed';

export interface TourSlide {
  title: string;
  description: string;
  icon: string;
  items: string[];
  goTo?: { panel: string; label: string };
}

export interface ChangelogEntry {
  version: string;
  date: string;
  items: { category: ChangeCategory; text: string }[];
  tourSlides?: TourSlide[];
}

export const CHANGELOG: ChangelogEntry[] = [
  {
    version: '1.1.0-dev',
    date: '2026-03-18',
    items: [
      { category: 'added', text: 'Python API — pip install flexfoil for scriptable airfoil analysis, polar sweeps, plotting, and a local web UI' },
      { category: 'added', text: 'Automatic L/D ratio column in Data Explorer, Plot Builder, and Polar plots' },
      { category: 'added', text: 'Custom computed columns — define algebraic expressions over existing data fields' },
      { category: 'added', text: 'Version & changelog dialogs accessible from the Help menu' },
      { category: 'added', text: 'Editable flap object model — add, edit, and remove multiple flaps with per-flap hinge y control' },
      { category: 'added', text: 'Reactive flap geometry — deflection changes update the airfoil instantly, no Apply button' },
      { category: 'added', text: 'Flap definitions persisted in run database (flaps_json column)' },
      { category: 'fixed', text: 'Flap hinge y-coordinate was computed at arc-length midpoint instead of the actual hinge x station' },
      { category: 'fixed', text: 'Flap geometry uses XFOIL-style fold trimming to eliminate surface self-intersection at the hinge' },
      { category: 'added', text: 'Multi-parameter sweep engine — sweep alpha, Reynolds, Mach, Ncrit, flap deflection, or flap hinge x/c' },
      { category: 'added', text: 'Matrix sweeps — combine two sweep parameters to generate a grid of polar series' },
      { category: 'added', text: 'Polar plot axes now include Reynolds, Mach, Ncrit, flap deflection, and flap hinge x/c' },
      { category: 'added', text: 'Solver queue in the status bar — see running sweeps with live progress bars and cancel from anywhere' },
      { category: 'added', text: 'Flap configurations included in shareable URLs' },
      { category: 'added', text: 'Solver result caching — repeated sweeps skip already-computed points and display cache hit stats' },
      { category: 'added', text: "What's New feature tour shown automatically on version update" },
      { category: 'added', text: 'Automatic IQR-based outlier removal toggle in Polar plot, Plot Builder, and Data Explorer' },
    ],
    tourSlides: [
      {
        title: 'Python API',
        description: 'pip install flexfoil — scriptable airfoil analysis from Python',
        icon: '🐍',
        items: [
          'Single-point solves and parallel polar sweeps from a Python script',
          'NACA generation, .dat loading, or custom coordinates',
          'Built-in Plotly and Matplotlib plotting',
          'flexfoil.serve() launches a local web UI sharing the same run database',
          'pandas DataFrames via polar.to_dataframe() and flexfoil.runs()',
        ],
      },
      {
        title: 'Flap Design System',
        description: 'Design multi-flap configurations with real-time geometry updates',
        icon: '✂',
        items: [
          'Add, edit, and remove multiple flaps with individual hinge control',
          'Deflection changes update the airfoil shape instantly',
          'XFOIL-style fold trimming for clean hinge surfaces',
          'Flap definitions saved with run data and included in shareable URLs',
        ],
        goTo: { panel: 'control', label: 'Open Geometry Control' },
      },
      {
        title: 'Multi-Parameter Sweeps',
        description: 'Sweep across multiple parameters in a single batch',
        icon: '⚡',
        items: [
          'Sweep alpha, Re, Mach, Ncrit, flap deflection, or hinge x/c',
          'Combine two parameters for matrix sweeps',
          'Polar plots support all swept parameters as axis choices',
          'Solver caching skips already-computed points across sweeps',
        ],
        goTo: { panel: 'solve', label: 'Open Solve panel' },
      },
      {
        title: 'Solver Status Bar',
        description: 'Monitor and manage all running solver jobs from the footer',
        icon: '🟢',
        items: [
          'Live status indicator — green when ready, pulsing yellow when running',
          'Click to open the queue popover with all recent jobs',
          'Per-job progress bars, elapsed time, and cancel buttons',
          'Works for polar sweeps, multi-sweeps, and inverse design jobs',
        ],
      },
      {
        title: 'Enhanced Data Analysis',
        description: 'New tools for exploring aerodynamic quantities',
        icon: '📐',
        items: [
          'Automatic L/D ratio in Data Explorer, Plot Builder, and Polar plots',
          'Custom computed columns with algebraic expressions over any data field',
        ],
        goTo: { panel: 'data-explorer', label: 'Open Data Explorer' },
      },
    ],
  },
  {
    version: '1.0.0',
    date: '2026-03-18',
    items: [
      { category: 'added', text: 'Initial public release of FlexFoil' },
      { category: 'added', text: 'Real-time airfoil analysis powered by RustFoil WASM solver' },
      { category: 'added', text: 'Airfoil library with NACA and Selig database' },
      { category: 'added', text: 'Interactive shape editing (camber/thickness splines)' },
      { category: 'added', text: 'Viscous and inviscid solver modes' },
      { category: 'added', text: 'Polar sweep with multi-series overlay' },
      { category: 'added', text: 'Flow visualization (streamlines, smoke, Cp, boundary layer)' },
      { category: 'added', text: 'Data Explorer with AG Grid and SPLOM correlogram' },
      { category: 'added', text: 'Plot Builder with scatter, line, bar, and histogram charts' },
      { category: 'added', text: 'SQLite run database with export/import' },
      { category: 'added', text: 'Mobile-responsive layout' },
    ],
  },
];
