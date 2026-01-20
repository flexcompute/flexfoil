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
  },
};

/**
 * Get a challenge by ID
 */
export function getChallenge(id: string): Challenge | undefined {
  return challenges[id];
}
