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
  showMeTourId?: string;
}

export interface ChangelogItem {
  category: ChangeCategory;
  text: string;
  showMeTourId?: string;
}

export interface ChangelogEntry {
  version: string;
  date: string;
  items: ChangelogItem[];
  tourSlides?: TourSlide[];
}

export const CHANGELOG: ChangelogEntry[] = [
  {
    version: '1.1.2',
    date: '2026-03-20',
    items: [
      { category: 'fixed', text: 'Tutorial no longer dismisses when clicking on panels or the overlay' },
      { category: 'changed', text: 'Hidden panel warning in tutorials replaced with a "Show Panel" button and animated cursor' },
    ],
  },
  {
    version: '1.1.1',
    date: '2026-03-19',
    items: [
      { category: 'fixed', text: 'File menu actions (New NACA, Import .dat, Export .dat, Export SVG) were disabled and not connected to anything' },
      { category: 'added', text: 'Export airfoil geometry as .dat (standard Selig format) or SVG from the File menu' },
      { category: 'added', text: 'Manually flag data points as outliers via right-click context menu in all plotting panels' },
      { category: 'added', text: 'Row grouping with aggregation in Data Explorer — group by airfoil, Re, Mach to see CLmax, CDmin, L/D_max per group' },
      { category: 'added', text: 'Custom aggregation functions: argmax/argmin for derived quantities like alpha_stall and alpha at L/D_max' },
      { category: 'added', text: 'Aerodynamic summary columns (α_stall, α @ L/D_max) that activate when row grouping is enabled' },
      { category: 'added', text: 'Aggregated data source in Plot Builder — plot group-level statistics (e.g. L/D_max vs Re) directly' },
      { category: 'added', text: 'Converged Only quick-filter button in Data Explorer to exclude non-converged points before aggregation' },
      { category: 'added', text: 'Python API: PolarResult now exposes cl_max, alpha_stall, ld_max, cd_min, and generic argmax/argmin/column_mean/column_median' },
      { category: 'added', text: 'Selig airfoil database browser with ~1,600 airfoils from the UIUC database — search and load any airfoil', showMeTourId: 'showMe:seligDatabase' },
      { category: 'added', text: 'Random Foil button to load a surprise airfoil from the Selig database' },
      { category: 'added', text: 'Smart Group button in Data Explorer — one click groups by polar configuration (airfoil + Re + Mach + Ncrit + flap) to see CL_max, L/D_max, α_stall per group', showMeTourId: 'showMe:smartGroup' },
      { category: 'added', text: 'Aggregated data source in Data Explorer correlogram and Polar Plot overlay — plot group-level statistics across all plotting surfaces' },
      { category: 'added', text: '"Show Me" interactive tutorials in the What\'s New dialog — guided walkthroughs for new features powered by the tour system' },
      { category: 'fixed', text: 'Data Explorer crashed with toFixed error when grouping by airfoil or other non-numeric columns' },
      { category: 'changed', text: 'Smart Group and correlogram data source now persist across view switches and page reloads' },
    ],
    tourSlides: [
      {
        title: 'Selig Airfoil Database',
        description: 'Browse ~1,600 real-world airfoils from the UIUC database',
        icon: '🔍',
        items: [
          'Search by name or designator — e387, clark, naca, and more',
          'Click any result to load it instantly onto the canvas',
          'Random Foil button picks a surprise airfoil for quick exploration',
          'Data sourced from the UIUC Airfoil Coordinates Database',
        ],
        goTo: { panel: 'library', label: 'Open Library' },
        showMeTourId: 'showMe:seligDatabase',
      },
      {
        title: 'Smart Group & Aggregation',
        description: 'One-click row grouping with aerodynamic summary statistics',
        icon: '📊',
        items: [
          'Smart Group button groups by airfoil + Re + Mach + Ncrit + flap configuration',
          'Summary columns: CL_max, CD_min, L/D_max, α_stall, α @ L/D_max',
          'Custom argmax/argmin aggregation for derived quantities',
          'Aggregated data source available in correlogram and Polar Plot overlay',
        ],
        goTo: { panel: 'data-explorer', label: 'Open Data Explorer' },
        showMeTourId: 'showMe:smartGroup',
      },
      {
        title: 'Manual Outlier Flagging',
        description: 'Right-click any data point to mark it as an outlier',
        icon: '🚩',
        items: [
          'Right-click context menu in all plotting panels (Polar, Plot Builder, Distributions)',
          'Flagged points are excluded from curve fits and aggregation',
          'Combine with the automatic IQR outlier filter for clean data',
          'Outlier flags persist in the run database',
        ],
      },
      {
        title: 'File Menu & Export',
        description: 'Export your airfoil geometry in standard formats',
        icon: '💾',
        items: [
          'Export as .dat (Selig format) for use in other tools',
          'Export as SVG for publication-quality figures',
          'New NACA and Import .dat actions now fully connected',
          'All actions accessible from the File menu',
        ],
      },
    ],
  },
  {
    version: '1.1.0-dev',
    date: '2026-03-18',
    items: [
      { category: 'added', text: 'Python API — pip install flexfoil for scriptable airfoil analysis, polar sweeps, plotting, and a local web UI' },
      { category: 'added', text: 'Automatic L/D ratio column in Data Explorer, Plot Builder, and Polar plots' },
      { category: 'added', text: 'Custom computed columns — define algebraic expressions over existing data fields', showMeTourId: 'showMe:dataAnalysis' },
      { category: 'added', text: 'Version & changelog dialogs accessible from the Help menu' },
      { category: 'added', text: 'Editable flap object model — add, edit, and remove multiple flaps with per-flap hinge y control', showMeTourId: 'showMe:flapDesign' },
      { category: 'added', text: 'Reactive flap geometry — deflection changes update the airfoil instantly, no Apply button' },
      { category: 'added', text: 'Flap definitions persisted in run database (flaps_json column)' },
      { category: 'fixed', text: 'Flap hinge y-coordinate was computed at arc-length midpoint instead of the actual hinge x station' },
      { category: 'fixed', text: 'Flap geometry uses XFOIL-style fold trimming to eliminate surface self-intersection at the hinge' },
      { category: 'added', text: 'Multi-parameter sweep engine — sweep alpha, Reynolds, Mach, Ncrit, flap deflection, or flap hinge x/c', showMeTourId: 'showMe:multiSweep' },
      { category: 'added', text: 'Matrix sweeps — combine two sweep parameters to generate a grid of polar series' },
      { category: 'added', text: 'Polar plot axes now include Reynolds, Mach, Ncrit, flap deflection, and flap hinge x/c' },
      { category: 'added', text: 'Solver queue in the status bar — see running sweeps with live progress bars and cancel from anywhere', showMeTourId: 'showMe:solverQueue' },
      { category: 'added', text: 'Flap configurations included in shareable URLs' },
      { category: 'added', text: 'Solver result caching — repeated sweeps skip already-computed points and display cache hit stats' },
      { category: 'added', text: "What's New feature tour shown automatically on version update" },
      { category: 'added', text: 'Automatic IQR-based outlier removal toggle in Polar plot, Plot Builder, and Data Explorer' },
      { category: 'added', text: 'Command palette (Cmd+K) — fuzzy search across panels, features, actions, and tutorials' },
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
        showMeTourId: 'showMe:flapDesign',
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
        showMeTourId: 'showMe:multiSweep',
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
        showMeTourId: 'showMe:solverQueue',
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
        showMeTourId: 'showMe:dataAnalysis',
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
