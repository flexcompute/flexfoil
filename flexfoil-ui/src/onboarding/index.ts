/**
 * Onboarding module - Barrel export
 */

export { OnboardingProvider, OnboardingContext } from './OnboardingProvider';
export type { OnboardingContextType } from './OnboardingProvider';
export { useOnboarding } from './useOnboarding';
export { tours, type TourId, type TourStep } from './tours';
export { challenges, getChallenge, type Challenge } from './challenges';
export * from './storage';
