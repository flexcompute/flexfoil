/**
 * useOnboarding - Hook for accessing onboarding context
 */

import { useContext } from 'react';
import { OnboardingContext, type OnboardingContextType } from './OnboardingProvider';

export function useOnboarding(): OnboardingContextType {
  const context = useContext(OnboardingContext);
  
  if (!context) {
    throw new Error('useOnboarding must be used within an OnboardingProvider');
  }
  
  return context;
}
