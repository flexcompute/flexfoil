/**
 * Challenge definitions for interactive tutorial steps
 * 
 * Challenges pause the tour until the user completes an action,
 * with an option to skip.
 */

import { useAirfoilStore } from '../stores/airfoilStore';
import { useVisualizationStore } from '../stores/visualizationStore';
import { useSolverJobStore } from '../stores/solverJobStore';
import { useCustomColumnStore } from '../stores/customColumnStore';

export interface Challenge {
  /** Unique identifier */
  id: string;
  /** Instruction shown to user */
  instruction: string;
  /** Hint text shown below instruction */
  hint?: string;
  /** Target value description (e.g., "10°") */
  targetValue?: string;
  /** Function that checks if challenge is complete */
  validate: () => boolean;
  /** Optional: highlight a specific element while challenge is active */
  highlightElement?: string;
  /** If true, the skip button is not shown - user must complete the challenge */
  noSkip?: boolean;
  /** Panel that must be visible for this challenge (will prompt to open if not visible) */
  requiredPanel?: string;
}

/**
 * Check if a panel is visible (exists and has non-zero dimensions)
 */
export function isPanelVisible(panelId: string): boolean {
  const panel = document.querySelector(`[data-tour="panel-${panelId}"]`);
  if (!panel) return false;
  const rect = panel.getBoundingClientRect();
  return rect.width > 0 && rect.height > 0;
}

/**
 * Check if any DOM element is visible (exists and has non-zero dimensions)
 */
export function isElementVisible(selector: string): boolean {
  const element = document.querySelector(selector);
  if (!element) return false;
  const rect = element.getBoundingClientRect();
  return rect.width > 0 && rect.height > 0;
}

/**
 * Get the display name for a panel
 */
export function getPanelDisplayName(panelId: string): string {
  const names: Record<string, string> = {
    'solve': 'Solve',
    'library': 'Library', 
    'visualization': 'Visualization',
    'control': 'Geometry Control',
    'control-mode': 'Geometry Control',
    'polar': 'Polar Plot',
    'canvas': 'Canvas',
    'spacing': 'Spacing',
    'properties': 'Properties',
    'info': 'Airfoil Info',
    'data-explorer': 'Data Explorer',
    'plot-builder': 'Plot Builder',
  };
  return names[panelId] || panelId;
}

/**
 * Map of element selectors to their parent panel IDs
 * Used to determine which panel to focus when an element isn't visible
 */
const ELEMENT_TO_PANEL: Record<string, string> = {
  // Control panel elements
  '[data-tour="thickness-slider"]': 'control',
  '[data-tour="camber-slider"]': 'control',
  '[data-tour="control-mode-parameters"]': 'control',
  '[data-tour="control-mode-camber"]': 'control',
  '[data-tour="control-mode-thickness"]': 'control',
  '[data-tour="panel-control"]': 'control',
  // Solve panel elements
  '[data-tour="solve-alpha"]': 'solve',
  '[data-tour="solve-run"]': 'solve',
  '[data-tour="solve-polar"]': 'solve',
  '[data-tour="panel-solve"]': 'solve',
  // Visualization panel elements
  '[data-tour="viz-streamlines"]': 'visualization',
  '[data-tour="viz-psi"]': 'visualization',
  '[data-tour="viz-smoke"]': 'visualization',
  '[data-tour="panel-visualization"]': 'visualization',
  // Geometry design elements
  '[data-tour="gdes-flaps"]': 'control',
  '[data-tour="gdes-add-flap"]': 'control',
  '[data-tour="control-mode-gdes"]': 'control',
  // Solver status
  '[data-tour="solver-status"]': 'solve',
  // Other panels
  '[data-tour="panel-library"]': 'library',
  '[data-tour="panel-spacing"]': 'spacing',
  '[data-tour="panel-polar"]': 'polar',
  '[data-tour="panel-canvas"]': 'canvas',
  '[data-tour="panel-properties"]': 'properties',
  // Data Explorer elements
  '[data-tour="de-tab-table"]': 'data-explorer',
  '[data-tour="de-tab-correlogram"]': 'data-explorer',
  '[data-tour="de-table-view"]': 'data-explorer',
  '[data-tour="de-export-csv"]': 'data-explorer',
  '[data-tour="de-export-db"]': 'data-explorer',
  '[data-tour="de-import-db"]': 'data-explorer',
  '[data-tour="de-column-chips"]': 'data-explorer',
  '[data-tour="de-encoding-controls"]': 'data-explorer',
  '[data-tour="de-correlogram-controls"]': 'data-explorer',
  '[data-tour="de-correlogram-plot"]': 'data-explorer',
  '[data-tour="panel-data-explorer"]': 'data-explorer',
  // Plot Builder elements
  '[data-tour="pb-tab-build"]': 'plot-builder',
  '[data-tour="pb-tab-saved"]': 'plot-builder',
  '[data-tour="pb-save-plot"]': 'plot-builder',
  '[data-tour="pb-controls"]': 'plot-builder',
  '[data-tour="pb-plot-area"]': 'plot-builder',
  '[data-tour="panel-plot-builder"]': 'plot-builder',
};

/**
 * Get the parent panel ID for a given element selector
 */
export function getParentPanel(selector: string): string | null {
  // Direct lookup
  if (ELEMENT_TO_PANEL[selector]) {
    return ELEMENT_TO_PANEL[selector];
  }
  
  // Try to find panel from element's data-tour attribute pattern
  // e.g., [data-tour="panel-solve"] -> solve
  const panelMatch = selector.match(/\[data-tour="panel-([^"]+)"\]/);
  if (panelMatch) {
    return panelMatch[1];
  }
  
  // Try to find by checking actual DOM element's parent
  const element = document.querySelector(selector);
  if (element) {
    // Walk up the DOM tree looking for a panel container
    let current = element.parentElement;
    while (current) {
      const tourAttr = current.getAttribute('data-tour');
      if (tourAttr?.startsWith('panel-')) {
        return tourAttr.replace('panel-', '');
      }
      current = current.parentElement;
    }
  }
  
  return null;
}

/**
 * Build HTML for element not visible warning with actionable "Show Panel" button
 * and an animated cursor to draw attention.
 */
export function buildElementNotVisibleHTML(panelId: string | null): string {
  const panelName = panelId ? getPanelDisplayName(panelId) : 'the required panel';
  const dataAttr = panelId ? ` data-open-panel="${panelId}"` : '';

  const cursorSvg = `<svg class="tour-animated-cursor" width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <path d="M5 3l14 8-6.5 1.5L11 19z" fill="currentColor" stroke="var(--bg-primary, #1a1a1a)" stroke-width="1.2" stroke-linejoin="round"/>
  </svg>`;

  return `
    <div class="tour-element-warning">
      <div class="tour-element-warning__header">
        ${cursorSvg}
        <span class="tour-element-warning__label">Panel hidden</span>
      </div>
      <div class="tour-element-warning__message">
        The <strong>${panelName}</strong> panel needs to be visible for this step.
      </div>
      <button class="tour-show-panel-btn"${dataAttr}>
        Show ${panelName}
      </button>
    </div>
  `;
}

/**
 * Get the current state for validation
 */
function getState() {
  const airfoil = useAirfoilStore.getState();
  const viz = useVisualizationStore.getState();
  const solverJobs = useSolverJobStore.getState();
  const customColumns = useCustomColumnStore.getState();
  return { airfoil, viz, solverJobs, customColumns };
}

/**
 * Pre-defined challenges that can be referenced by ID in tour steps
 */
export const challenges: Record<string, Challenge> = {
  // Solver challenges
  'set-alpha-10': {
    id: 'set-alpha-10',
    instruction: 'Set the angle of attack to 10°',
    hint: 'Type 10 in the Alpha input field in the Solve panel',
    targetValue: '10°',
    validate: () => {
      const { airfoil } = getState();
      return Math.abs(airfoil.displayAlpha - 10) < 0.5;
    },
    highlightElement: '[data-tour="solve-alpha"]',
    requiredPanel: 'solve',
  },

  'run-analysis': {
    id: 'run-analysis',
    instruction: 'Click the Run button to analyze',
    hint: 'This computes lift and moment coefficients',
    validate: () => {
      // For now, we'll skip validation since analysis results aren't stored in global state
      // The user just needs to click the button - we trust they did it
      // In the future, we could add a "lastAnalysisTimestamp" to the store
      return true; // Auto-complete this challenge
    },
    highlightElement: '[data-tour="solve-run"]',
    requiredPanel: 'solve',
  },

  'generate-polar': {
    id: 'generate-polar',
    instruction: 'Generate a polar sweep',
    hint: 'Click "Generate Polar" to compute CL across multiple angles',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.polarData.some(s => s.points.length > 5);
    },
    highlightElement: '[data-tour="solve-polar"]',
    requiredPanel: 'solve',
  },

  // Library challenges
  'load-naca-2412': {
    id: 'load-naca-2412',
    instruction: 'Load the NACA 2412 airfoil',
    hint: 'Click the 2412 preset button or type 2412 and click Generate',
    targetValue: 'NACA 2412',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.name.includes('2412');
    },
    highlightElement: '[data-tour="panel-library"]',
    requiredPanel: 'library',
  },

  'change-airfoil': {
    id: 'change-airfoil',
    instruction: 'Load any different airfoil',
    hint: 'Try a preset like 4412 or generate your own NACA code',
    validate: () => {
      const { airfoil } = getState();
      // Check if name changed from default
      return !airfoil.name.includes('0012');
    },
    highlightElement: '[data-tour="panel-library"]',
    requiredPanel: 'library',
  },

  // Parameter slider challenges
  'adjust-thickness': {
    id: 'adjust-thickness',
    instruction: 'Adjust the Thickness Scale slider',
    hint: 'Drag the slider to make the airfoil thicker or thinner',
    validate: () => {
      const { airfoil } = getState();
      // Check if thickness scale has been changed from default (1.0)
      return Math.abs(airfoil.thicknessScale - 1.0) > 0.05;
    },
    highlightElement: '[data-tour="thickness-slider"]',
    requiredPanel: 'control',
  },

  'adjust-camber': {
    id: 'adjust-camber',
    instruction: 'Adjust the Camber Scale slider',
    hint: 'Drag the slider to change the airfoil curvature',
    validate: () => {
      const { airfoil } = getState();
      // Check if camber scale has been changed from default (1.0)
      return Math.abs(airfoil.camberScale - 1.0) > 0.05;
    },
    highlightElement: '[data-tour="camber-slider"]',
    requiredPanel: 'control',
  },

  // Spline control mode challenges (for advanced editing tour)
  'switch-to-camber': {
    id: 'switch-to-camber',
    instruction: 'Switch to Camber editing mode',
    hint: 'Click the "Camber" button in Geometry Control panel',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.controlMode === 'camber-spline';
    },
    highlightElement: '[data-tour="control-mode-camber"]',
    requiredPanel: 'control-mode',
  },

  'modify-camber': {
    id: 'modify-camber',
    instruction: 'Drag a camber control point on the canvas',
    hint: 'Click and drag any blue control point up or down to change the camber line',
    validate: () => {
      const { airfoil } = getState();
      // Check if in camber mode and any middle control point has been moved
      if (airfoil.controlMode !== 'camber-spline') return false;
      // Check if any control point (not at endpoints) has non-zero camber
      const middlePoints = airfoil.camberControlPoints.filter(
        (p) => p.x > 0.05 && p.x < 0.95
      );
      // If any middle point has |y| > 0.005, user has modified the camber
      return middlePoints.some((p) => Math.abs(p.y) > 0.005);
    },
    highlightElement: '[data-tour="panel-canvas"]',
    requiredPanel: 'control',
  },

  'switch-to-thickness': {
    id: 'switch-to-thickness',
    instruction: 'Switch to Thickness editing mode',
    hint: 'Click the "Thickness" button in Geometry Control panel',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.controlMode === 'thickness-spline';
    },
    highlightElement: '[data-tour="control-mode-thickness"]',
    requiredPanel: 'control-mode',
  },

  'modify-thickness': {
    id: 'modify-thickness',
    instruction: 'Drag a thickness control point on the canvas',
    hint: 'Click and drag any orange control point up or down to change thickness',
    validate: () => {
      const { airfoil } = getState();
      // Check if in thickness mode
      if (airfoil.controlMode !== 'thickness-spline') return false;
      // Check if the max thickness differs from the original ~0.06 (NACA 0012 is 12% thick, half-thickness ~0.06)
      const maxT = Math.max(...airfoil.thicknessControlPoints.map((p) => p.t));
      // Original NACA 0012 max half-thickness is about 0.06, allow some tolerance
      return Math.abs(maxT - 0.06) > 0.005;
    },
    highlightElement: '[data-tour="panel-canvas"]',
    requiredPanel: 'control',
  },

  // Panel visibility challenges
  'open-visualization-panel': {
    id: 'open-visualization-panel',
    instruction: 'Open the Visualization panel',
    hint: 'Go to Window menu and click on "Visualization", or click its tab if visible',
    validate: () => {
      // Check if the visualization panel element exists AND is visible (active tab)
      const vizPanel = document.querySelector('[data-tour="panel-visualization"]');
      if (!vizPanel) return false;
      
      // Check if the element is actually visible (not hidden by tab)
      // FlexLayout hides inactive tabs with display:none on the parent
      const rect = vizPanel.getBoundingClientRect();
      return rect.width > 0 && rect.height > 0;
    },
    highlightElement: '[data-tour="menu-window"]',
  },

  'open-solve-panel': {
    id: 'open-solve-panel',
    instruction: 'Open the Solve panel',
    hint: 'Go to Window menu and click on "Solve", or click its tab if visible',
    validate: () => {
      // Check if the solve panel element exists AND is visible (active tab)
      const solvePanel = document.querySelector('[data-tour="panel-solve"]');
      if (!solvePanel) return false;
      
      // Check if the element is actually visible (not hidden by tab)
      const rect = solvePanel.getBoundingClientRect();
      return rect.width > 0 && rect.height > 0;
    },
    highlightElement: '[data-tour="menu-window"]',
  },

  // Visualization challenges
  'enable-streamlines': {
    id: 'enable-streamlines',
    instruction: 'Enable streamline visualization',
    hint: 'Toggle streamlines on in the Visualization panel',
    validate: () => {
      const { viz } = getState();
      return viz.showStreamlines;
    },
    highlightElement: '[data-tour="viz-streamlines"]',
    requiredPanel: 'visualization',
  },

  'enable-psi': {
    id: 'enable-psi',
    instruction: 'Enable stream function visualization',
    hint: 'Toggle Stream Function (ψ) on',
    validate: () => {
      const { viz } = getState();
      return viz.showPsiContours;
    },
    highlightElement: '[data-tour="viz-psi"]',
    requiredPanel: 'visualization',
  },

  'enable-smoke': {
    id: 'enable-smoke',
    instruction: 'Enable smoke visualization',
    hint: 'Toggle Smoke on to see animated flow',
    validate: () => {
      const { viz } = getState();
      return viz.showSmoke;
    },
    highlightElement: '[data-tour="viz-smoke"]',
    requiredPanel: 'visualization',
  },

  // --- Flap challenges ---

  'add-flap': {
    id: 'add-flap',
    instruction: 'Add a flap to the airfoil',
    hint: 'Click the "Add Flap" button in the Geometry tab',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.geometryDesign.flaps.length > 0;
    },
    highlightElement: '[data-tour="gdes-add-flap"]',
    requiredPanel: 'control',
  },

  'deflect-flap': {
    id: 'deflect-flap',
    instruction: 'Deflect the flap',
    hint: 'Drag the deflection slider or type a value (e.g. 10)',
    targetValue: 'any non-zero deflection',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.geometryDesign.flaps.some(f => Math.abs(f.deflection) > 0.5);
    },
    highlightElement: '[data-tour="gdes-flaps"]',
    requiredPanel: 'control',
  },

  // --- Sweep challenges ---

  'configure-matrix-sweep': {
    id: 'configure-matrix-sweep',
    instruction: 'Enable a secondary sweep parameter',
    hint: 'Click "Add Secondary Sweep" in the polar section of the Solve panel',
    validate: () => {
      const rows = document.querySelectorAll('[data-tour="solve-polar"] .sweep-param-row');
      return rows.length >= 2;
    },
    highlightElement: '[data-tour="solve-polar"]',
    requiredPanel: 'solve',
  },

  'run-sweep': {
    id: 'run-sweep',
    instruction: 'Run a polar sweep',
    hint: 'Click "Generate Polar" to start the sweep',
    validate: () => {
      const { solverJobs } = getState();
      return solverJobs.jobs.length > 0;
    },
    highlightElement: '[data-tour="solve-polar"]',
    requiredPanel: 'solve',
  },

  // --- Data analysis challenges ---

  'add-computed-column': {
    id: 'add-computed-column',
    instruction: 'Create a custom computed column',
    hint: 'Click "Add Column" in the Data Explorer and enter an expression like cl/cd',
    validate: () => {
      const { customColumns } = getState();
      return customColumns.columns.length > 0;
    },
    highlightElement: '[data-tour="de-column-chips"]',
    requiredPanel: 'data-explorer',
  },

  // --- Control mode challenges ---

  'switch-to-gdes': {
    id: 'switch-to-gdes',
    instruction: 'Switch to Geometry editing mode',
    hint: 'Click the "Geometry" button in Geometry Control panel',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.controlMode === 'geometry-design';
    },
    highlightElement: '[data-tour="control-mode-gdes"]',
    requiredPanel: 'control',
  },
};

/**
 * Get a challenge by ID
 */
export function getChallenge(id: string): Challenge | undefined {
  return challenges[id];
}
