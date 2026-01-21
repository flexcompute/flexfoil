/**
 * Onboarding state persistence using localStorage
 */

const STORAGE_KEY = 'flexfoil-onboarding';
const CURRENT_VERSION = '1.0.0';

export interface OnboardingState {
  completedTours: string[];
  /** Progress within each tour (step index) - allows resuming */
  tourProgress: Record<string, number>;
  lastSeenVersion: string;
}

const defaultState: OnboardingState = {
  completedTours: [],
  tourProgress: {},
  lastSeenVersion: CURRENT_VERSION,
};

export function getOnboardingState(): OnboardingState {
  try {
    const stored = localStorage.getItem(STORAGE_KEY);
    if (!stored) return defaultState;
    
    const state = JSON.parse(stored) as OnboardingState;
    
    // If version changed, reset tours so users see new content
    if (state.lastSeenVersion !== CURRENT_VERSION) {
      return { ...defaultState, lastSeenVersion: CURRENT_VERSION };
    }
    
    // Ensure tourProgress exists (migration from old format)
    if (!state.tourProgress) {
      state.tourProgress = {};
    }
    
    return state;
  } catch {
    return defaultState;
  }
}

export function markTourComplete(tourId: string): void {
  const state = getOnboardingState();
  if (!state.completedTours.includes(tourId)) {
    state.completedTours.push(tourId);
  }
  // Clear progress when tour is completed
  delete state.tourProgress[tourId];
  localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
}

export function hasCompletedTour(tourId: string): boolean {
  const state = getOnboardingState();
  return state.completedTours.includes(tourId);
}

/** Save the current step index for a tour (for resuming later) */
export function saveTourProgress(tourId: string, stepIndex: number): void {
  const state = getOnboardingState();
  state.tourProgress[tourId] = stepIndex;
  localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
}

/** Get the saved step index for a tour (returns 0 if none saved) */
export function getTourProgress(tourId: string): number {
  const state = getOnboardingState();
  return state.tourProgress[tourId] ?? 0;
}

/** Check if a tour has any saved progress (even at step 0) */
export function hasTourProgress(tourId: string): boolean {
  const state = getOnboardingState();
  return tourId in state.tourProgress;
}

/** Clear progress for a specific tour (when restarting from beginning) */
export function clearTourProgress(tourId: string): void {
  const state = getOnboardingState();
  delete state.tourProgress[tourId];
  localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
}

export function resetOnboarding(): void {
  localStorage.removeItem(STORAGE_KEY);
}
