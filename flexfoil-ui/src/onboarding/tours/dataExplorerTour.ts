/**
 * Data Explorer & Correlogram Tour
 *
 * Walks users through the Data Explorer table, correlogram,
 * encoding controls, right-click context menu, and saving
 * plots via the Plot Builder.
 *
 * Triggered from the ? button in the Data Explorer header
 * or from the Help menu.
 */

import type { TourStep } from './index';
import { useRouteUiStore } from '../../stores/routeUiStore';

const tipBox = (text: string) => `
  <div style="margin-top: 12px; padding: 10px; background: var(--bg-tertiary); border-radius: 6px; border-left: 3px solid var(--accent-secondary);">
    <strong>Tip:</strong> ${text}
  </div>`;

export const dataExplorerTour: TourStep[] = [
  // ── Intro ──
  {
    element: '[data-tour="panel-data-explorer"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Data Explorer',
      description: `This panel is your home for browsing, filtering, and visualising every solver run you've made. It has two views: a <strong>Table</strong> and a <strong>Correlogram</strong>.${tipBox('You can reopen this guide anytime from the <strong>Help</strong> menu.')}`,
      side: 'left',
      align: 'start',
    },
  },

  // ── Table view ──
  {
    element: '[data-tour="de-tab-table"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Table View',
      description: 'The <strong>Table</strong> tab shows every run as a row in a full-featured data grid. You can sort, filter, group, and even create AG Grid charts directly from here.',
      side: 'bottom',
      align: 'start',
    },
    onBeforeHighlight: () => {
      useRouteUiStore.getState().setDataExplorerView('table');
    },
  },
  {
    element: '[data-tour="de-table-view"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'AG Grid Features',
      description: `The table supports enterprise-grade features:
        <ul style="margin: 8px 0 0 0; padding-left: 18px; font-size: 12px; line-height: 1.6;">
          <li><strong>Column filters</strong> — click the ≡ icon on any header</li>
          <li><strong>Row grouping</strong> — drag a column header to the "Row Group" bar at the top</li>
          <li><strong>Column sidebar</strong> — toggle columns on/off from the right sidebar</li>
          <li><strong>Cell selection</strong> — click and drag to select cells for aggregation in the status bar</li>
        </ul>`,
      side: 'left',
      align: 'center',
    },
  },
  {
    element: '[data-tour="de-export-csv"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Export CSV',
      description: 'Download the currently visible (filtered) rows as a <strong>.csv</strong> file you can open in Excel, Google Sheets, or any data tool.',
      side: 'bottom',
      align: 'end',
    },
  },
  {
    element: '[data-tour="de-export-db"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Export / Import Database',
      description: 'Export the entire run database as a <strong>.sqlite</strong> file for backup, or import one to restore a previous session. Use <strong>Import DB</strong> to load a previously exported database.',
      side: 'bottom',
      align: 'end',
    },
  },

  // ── Switch to Correlogram ──
  {
    element: '[data-tour="de-tab-correlogram"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Correlogram View',
      description: 'Switch to the <strong>Correlogram</strong> to see a scatter-plot matrix of all your numeric fields. Each off-diagonal cell plots one variable against another; the diagonal shows histograms.',
      side: 'bottom',
      align: 'start',
    },
    onBeforeHighlight: () => {
      useRouteUiStore.getState().setDataExplorerView('correlogram');
    },
  },

  // ── Column chips ──
  {
    element: '[data-tour="de-column-chips"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Choose Columns',
      description: 'Toggle which numeric fields appear in the matrix. Each active chip becomes a row and column. You need at least two.',
      side: 'bottom',
      align: 'start',
    },
  },

  // ── Encoding controls ──
  {
    element: '[data-tour="de-encoding-controls"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Encoding Controls',
      description: `Map visual properties to your data:
        <ul style="margin: 8px 0 0 0; padding-left: 18px; font-size: 12px; line-height: 1.6;">
          <li><strong>Group By</strong> — split traces by a categorical field (or use Auto)</li>
          <li><strong>Color</strong> — color points by a numeric or categorical field</li>
          <li><strong>Marker Size / Type</strong> — encode additional dimensions</li>
        </ul>
        ${tipBox('The <strong>Auto (Smart)</strong> grouping detects sweep series automatically.')}`,
      side: 'bottom',
      align: 'center',
    },
  },

  // ── The SPLOM itself ──
  {
    element: '[data-tour="de-correlogram-plot"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Reading the Correlogram',
      description: `<strong>Diagonal</strong> cells are histograms showing each variable's distribution. <strong>Off-diagonal</strong> cells are scatter plots, annotated with the Pearson <em>r</em> correlation coefficient.
        ${tipBox('<strong>Click any scatter point</strong> to restore that run\'s historical flowfield on the canvas.')}`,
      side: 'left',
      align: 'center',
    },
  },

  // ── Right-click context menu ──
  {
    element: '[data-tour="de-correlogram-plot"]',
    focusPanel: 'data-explorer',
    popover: {
      title: 'Open in Plot Builder',
      description: `<strong>Right-click</strong> any off-diagonal scatter cell to open a context menu with <strong>"Open in Plot Builder"</strong>. This pre-fills the X/Y axes and switches you to the Plot Builder for a full-size, single-chart view.
        <div style="margin-top: 10px; font-size: 12px; color: var(--text-secondary);">
          This is the quickest way to go from exploratory matrix → focused plot.
        </div>`,
      side: 'left',
      align: 'center',
    },
  },

  // ── Transition to Plot Builder ──
  {
    element: '[data-tour="panel-plot-builder"]',
    focusPanel: 'plot-builder',
    popover: {
      title: 'Plot Builder',
      description: 'The Plot Builder gives you a dedicated space to create single charts from your run data. Let\'s see how it works.',
      side: 'left',
      align: 'start',
    },
  },
  {
    element: '[data-tour="pb-controls"]',
    focusPanel: 'plot-builder',
    popover: {
      title: 'Build Controls',
      description: `Pick your chart type (scatter, line, bar, histogram), set X/Y axes with linear or log scale, choose data source, grouping, and visual encodings — all in one toolbar.`,
      side: 'bottom',
      align: 'start',
    },
  },
  {
    element: '[data-tour="pb-plot-area"]',
    focusPanel: 'plot-builder',
    popover: {
      title: 'Interactive Plot',
      description: 'Your chart renders here. Use the Plotly toolbar to zoom, pan, or download as PNG. <strong>Click a point</strong> to restore that run\'s flowfield, just like in the correlogram.',
      side: 'left',
      align: 'center',
    },
  },

  // ── Saving plots ──
  {
    element: '[data-tour="pb-save-plot"]',
    focusPanel: 'plot-builder',
    popover: {
      title: 'Save a Plot',
      description: 'Happy with your chart? Click <strong>Save Plot</strong> to store the current configuration (axes, chart type, encodings) for later. You\'ll be prompted for a title — or just press Enter for the default.',
      side: 'bottom',
      align: 'end',
    },
  },
  {
    element: '[data-tour="pb-tab-saved"]',
    focusPanel: 'plot-builder',
    popover: {
      title: 'Saved Plots',
      description: 'Switch to the <strong>Saved</strong> tab to browse all your saved plot configurations. Click any card to reload it instantly, or use the ✕ to delete.',
      side: 'bottom',
      align: 'start',
    },
  },

  // ── Wrap-up ──
  {
    popover: {
      title: 'That\'s the Data Explorer!',
      description: `Here's a quick recap:
        <ul style="margin: 8px 0 0 0; padding-left: 18px; font-size: 12px; line-height: 1.6;">
          <li><strong>Table</strong> — filter, sort, group, and export runs</li>
          <li><strong>Correlogram</strong> — spot trends across variables at a glance</li>
          <li><strong>Right-click → Plot Builder</strong> — zoom in on any relationship</li>
          <li><strong>Save Plot</strong> — bookmark configurations for later</li>
          <li><strong>Click a point</strong> — restore any historical flowfield</li>
        </ul>`,
      side: 'over',
      align: 'center',
    },
  },
];
