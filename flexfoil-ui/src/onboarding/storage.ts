/**
 * Onboarding state persistence using localStorage
 */

const STORAGE_KEY = 'flexfoil-onboarding';
const CURRENT_VERSION = '1.0.0';

export interface OnboardingState {
  completedTours: string[];
  lastSeenVersion: string;
}

const defaultState: OnboardingState = {
  completedTours: [],
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
  localStorage.setItem(STORAGE_KEY, JSON.stringify(state));
}

export function hasCompletedTour(tourId: string): boolean {
  const state = getOnboardingState();
  return state.completedTours.includes(tourId);
}

export function resetOnboarding(): void {
  localStorage.removeItem(STORAGE_KEY);
}
