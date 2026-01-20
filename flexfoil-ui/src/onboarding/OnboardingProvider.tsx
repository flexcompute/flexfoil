/**
 * OnboardingProvider - Context provider for tour management
 * 
 * Wraps the app and provides tour control functions via context.
 * Uses driver.js for the actual tour UI.
 * 
 * Supports interactive challenges where users must complete an action
 * before proceeding (with skip option).
 */

import { createContext, useCallback, useRef, useState, useEffect, type ReactNode } from 'react';
import { driver, type Driver, type DriveStep } from 'driver.js';
import 'driver.js/dist/driver.css';

import { tours, type TourId, type TourStep } from './tours';
import { getChallenge, type Challenge } from './challenges';
import { markTourComplete, hasCompletedTour as checkTourCompleted, resetOnboarding } from './storage';

export interface OnboardingContextType {
  /** Start a specific tour */
  startTour: (tourId: TourId) => void;
  /** End the current tour */
  endTour: () => void;
  /** Check if a tour has been completed */
  hasCompletedTour: (tourId: TourId) => boolean;
  /** Reset all tour progress */
  resetAllTours: () => void;
  /** Whether a tour is currently active */
  isActive: boolean;
  /** The currently running tour ID */
  currentTour: TourId | null;
}

export const OnboardingContext = createContext<OnboardingContextType | null>(null);

interface OnboardingProviderProps {
  children: ReactNode;
}

export function OnboardingProvider({ children }: OnboardingProviderProps) {
  const [isActive, setIsActive] = useState(false);
  const [currentTour, setCurrentTour] = useState<TourId | null>(null);
  const driverRef = useRef<Driver | null>(null);
  const challengeIntervalRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const currentChallengeRef = useRef<Challenge | null>(null);

  // Clean up challenge polling
  const clearChallengePolling = useCallback(() => {
    if (challengeIntervalRef.current) {
      clearInterval(challengeIntervalRef.current);
      challengeIntervalRef.current = null;
    }
    currentChallengeRef.current = null;
  }, []);

  // Build challenge HTML for popover
  const buildChallengeHTML = useCallback((challenge: Challenge, isComplete: boolean) => {
    return `
      <div class="tour-challenge ${isComplete ? 'tour-challenge--complete' : ''}">
        <div class="tour-challenge__header">
          <span class="tour-challenge__icon">${isComplete ? '✓' : '→'}</span>
          <span class="tour-challenge__label">${isComplete ? 'Complete!' : 'Try it'}</span>
        </div>
        <div class="tour-challenge__instruction">${challenge.instruction}</div>
        ${challenge.targetValue ? `<div class="tour-challenge__target">Target: <strong>${challenge.targetValue}</strong></div>` : ''}
        ${challenge.hint && !isComplete ? `<div class="tour-challenge__hint">${challenge.hint}</div>` : ''}
        ${isComplete ? '<div class="tour-challenge__success">Great job! Click Next to continue.</div>' : ''}
      </div>
    `;
  }, []);

  // Convert TourStep to DriveStep, injecting challenge UI
  const convertToDriverSteps = useCallback((tourSteps: TourStep[]): DriveStep[] => {
    return tourSteps.map((step) => {
      const challenge = step.challengeId ? getChallenge(step.challengeId) : null;
      
      // Build description with challenge UI if needed
      let description = step.popover?.description ?? '';
      if (challenge) {
        const isComplete = challenge.validate();
        description = `${description}${buildChallengeHTML(challenge, isComplete)}`;
      }

      return {
        element: step.element,
        popover: {
          ...step.popover,
          description,
        },
      };
    });
  }, [buildChallengeHTML]);

  // Initialize driver instance with theme-aware styling
  const getDriverInstance = useCallback((tourSteps: TourStep[]) => {
    if (driverRef.current) {
      driverRef.current.destroy();
    }
    clearChallengePolling();

    // Check theme for overlay color
    const isDark = document.documentElement.getAttribute('data-theme') !== 'light';

    driverRef.current = driver({
      showProgress: true,
      showButtons: ['next', 'previous', 'close'],
      nextBtnText: 'Next',
      prevBtnText: 'Back',
      doneBtnText: 'Done',
      progressText: '{{current}} of {{total}}',
      allowClose: true,
      overlayColor: isDark ? 'rgba(0, 0, 0, 0.75)' : 'rgba(0, 0, 0, 0.5)',
      stagePadding: 8,
      stageRadius: 8,
      popoverClass: 'flexfoil-tour-popover',
      
      onHighlightStarted: () => {
        clearChallengePolling();
        
        // Find the original tour step to check for challenge
        const stepIndex = driverRef.current?.getActiveIndex() ?? 0;
        const originalStep = tourSteps[stepIndex];
        
        if (originalStep?.challengeId) {
          const challenge = getChallenge(originalStep.challengeId);
          if (challenge) {
            currentChallengeRef.current = challenge;
            
            // Start polling for challenge completion
            challengeIntervalRef.current = setInterval(() => {
              if (challenge.validate()) {
                clearChallengePolling();
                
                // Update the popover to show completion
                const popoverEl = document.querySelector('.driver-popover-description');
                if (popoverEl && originalStep.popover) {
                  const baseDesc = originalStep.popover.description ?? '';
                  popoverEl.innerHTML = `${baseDesc}${buildChallengeHTML(challenge, true)}`;
                }
                
                // Enable the next button
                const nextBtn = document.querySelector('.driver-popover-next-btn') as HTMLButtonElement;
                if (nextBtn) {
                  nextBtn.disabled = false;
                  nextBtn.classList.remove('driver-popover-btn--disabled');
                }
              }
            }, 200);
          }
        }
      },

      onPopoverRender: (popover) => {
        // Check if current step has an incomplete challenge
        const stepIndex = driverRef.current?.getActiveIndex() ?? 0;
        const originalStep = tourSteps[stepIndex];
        
        if (originalStep?.challengeId) {
          const challenge = getChallenge(originalStep.challengeId);
          if (challenge && !challenge.validate()) {
            // Add skip button and disable next
            const footer = popover.footerButtons;
            const nextBtn = footer.querySelector('.driver-popover-next-btn') as HTMLButtonElement;
            
            if (nextBtn) {
              // Create skip button
              const skipBtn = document.createElement('button');
              skipBtn.className = 'driver-popover-skip-btn';
              skipBtn.textContent = 'Skip';
              skipBtn.onclick = () => {
                clearChallengePolling();
                driverRef.current?.moveNext();
              };
              
              // Insert skip before next
              nextBtn.parentNode?.insertBefore(skipBtn, nextBtn);
              
              // Disable next until challenge complete
              nextBtn.disabled = true;
              nextBtn.classList.add('driver-popover-btn--disabled');
            }
          }
        }
      },

      onDestroyed: () => {
        clearChallengePolling();
        if (currentTour) {
          markTourComplete(currentTour);
        }
        setIsActive(false);
        setCurrentTour(null);
      },

      onCloseClick: () => {
        clearChallengePolling();
        if (currentTour) {
          markTourComplete(currentTour);
        }
        driverRef.current?.destroy();
      },
    });

    return driverRef.current;
  }, [currentTour, clearChallengePolling, buildChallengeHTML]);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      clearChallengePolling();
      if (driverRef.current) {
        driverRef.current.destroy();
      }
    };
  }, [clearChallengePolling]);

  const startTour = useCallback((tourId: TourId) => {
    const tourSteps = tours[tourId];
    if (!tourSteps) {
      console.warn(`Tour "${tourId}" not found`);
      return;
    }

    setCurrentTour(tourId);
    setIsActive(true);

    // Small delay to ensure DOM is ready
    setTimeout(() => {
      const driverInstance = getDriverInstance(tourSteps);
      const driverSteps = convertToDriverSteps(tourSteps);
      driverInstance.setSteps(driverSteps);
      driverInstance.drive();
    }, 100);
  }, [getDriverInstance, convertToDriverSteps]);

  const endTour = useCallback(() => {
    clearChallengePolling();
    if (driverRef.current) {
      if (currentTour) {
        markTourComplete(currentTour);
      }
      driverRef.current.destroy();
    }
    setIsActive(false);
    setCurrentTour(null);
  }, [currentTour, clearChallengePolling]);

  const hasCompletedTour = useCallback((tourId: TourId): boolean => {
    return checkTourCompleted(tourId);
  }, []);

  const resetAllTours = useCallback(() => {
    resetOnboarding();
  }, []);

  const contextValue: OnboardingContextType = {
    startTour,
    endTour,
    hasCompletedTour,
    resetAllTours,
    isActive,
    currentTour,
  };

  return (
    <OnboardingContext.Provider value={contextValue}>
      {children}
    </OnboardingContext.Provider>
  );
}
