/**
 * Challenge definitions for interactive tutorial steps
 * 
 * Challenges pause the tour until the user completes an action,
 * with an option to skip.
 */

import { useAirfoilStore } from '../stores/airfoilStore';
import { useVisualizationStore } from '../stores/visualizationStore';

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
 * Get the display name for a panel
 */
export function getPanelDisplayName(panelId: string): string {
  const names: Record<string, string> = {
    'solve': 'Solve',
    'library': 'Library', 
    'visualization': 'Visualization',
    'control-mode': 'Control Mode',
    'polar': 'Polar Plot',
    'info': 'Airfoil Info',
  };
  return names[panelId] || panelId;
}

/**
 * Get the current state for validation
 */
function getState() {
  const airfoil = useAirfoilStore.getState();
  const viz = useVisualizationStore.getState();
  return { airfoil, viz };
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
      return airfoil.polarData.length > 5;
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

  // Control mode challenges
  'switch-to-camber': {
    id: 'switch-to-camber',
    instruction: 'Switch to Camber editing mode',
    hint: 'Click the "Camber" button in Control Mode panel',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.controlMode === 'camber-spline';
    },
    highlightElement: '[data-tour="control-mode-camber"]',
    requiredPanel: 'control-mode',
  },

  'switch-to-thickness': {
    id: 'switch-to-thickness',
    instruction: 'Switch to Thickness editing mode',
    hint: 'Click the "Thickness" button in Control Mode panel',
    validate: () => {
      const { airfoil } = getState();
      return airfoil.controlMode === 'thickness-spline';
    },
    highlightElement: '[data-tour="control-mode-thickness"]',
    requiredPanel: 'control-mode',
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
    noSkip: true,  // Must complete - subsequent steps depend on this panel being visible
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
    noSkip: true,  // Must complete - subsequent steps depend on this panel being visible
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
};

/**
 * Get a challenge by ID
 */
export function getChallenge(id: string): Challenge | undefined {
  return challenges[id];
}
