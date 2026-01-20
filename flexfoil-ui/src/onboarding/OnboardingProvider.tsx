/**
 * OnboardingProvider - Context provider for tour management
 * 
 * Wraps the app and provides tour control functions via context.
 * Uses driver.js for the actual tour UI.
 * 
 * Supports interactive challenges where users must complete an action
 * before proceeding (with skip option).
 * 
 * Supports focusing panels before highlighting steps.
 */

import { createContext, useCallback, useRef, useState, useEffect, type ReactNode } from 'react';
import { driver, type Driver, type DriveStep } from 'driver.js';
import 'driver.js/dist/driver.css';

import { tours, type TourId, type TourStep } from './tours';
import { getChallenge, isPanelVisible, getPanelDisplayName, type Challenge } from './challenges';
import { markTourComplete, hasCompletedTour as checkTourCompleted, resetOnboarding, saveTourProgress, getTourProgress, clearTourProgress } from './storage';
import { useLayout } from '../contexts/LayoutContext';

export interface OnboardingContextType {
  /** Start a specific tour (resumes from last position unless fromBeginning=true) */
  startTour: (tourId: TourId, fromBeginning?: boolean) => void;
  /** End the current tour */
  endTour: () => void;
  /** Check if a tour has been completed (all steps finished) */
  hasCompletedTour: (tourId: TourId) => boolean;
  /** Check if a tour has been started (has progress or is complete) - use for auto-start checks */
  hasStartedTour: (tourId: TourId) => boolean;
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
  
  // Get layout context for panel focusing
  const { openPanel } = useLayout();

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
    // Check if required panel is visible
    const panelVisible = challenge.requiredPanel ? isPanelVisible(challenge.requiredPanel) : true;
    const panelName = challenge.requiredPanel ? getPanelDisplayName(challenge.requiredPanel) : '';
    
    // If panel not visible, show panel requirement first
    const panelWarning = !panelVisible && challenge.requiredPanel ? `
      <div class="tour-challenge__panel-warning">
        <span class="tour-challenge__icon">⚠</span>
        <span>First, open the <strong>${panelName}</strong> panel from the Window menu</span>
      </div>
    ` : '';
    
    return `
      <div class="tour-challenge ${isComplete ? 'tour-challenge--complete' : ''}">
        <div class="tour-challenge__header">
          <span class="tour-challenge__icon">${isComplete ? '✓' : '→'}</span>
          <span class="tour-challenge__label">${isComplete ? 'Complete!' : 'Try it'}</span>
        </div>
        ${panelWarning}
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
  const getDriverInstance = useCallback((tourSteps: TourStep[], tourId: TourId) => {
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
      allowClose: false,  // Don't close on overlay click - use X button only
      overlayClickNext: false,  // Don't advance on overlay click either
      overlayColor: isDark ? 'rgba(0, 0, 0, 0.75)' : 'rgba(0, 0, 0, 0.5)',
      stagePadding: 8,
      stageRadius: 8,
      popoverClass: 'flexfoil-tour-popover',
      
      onHighlightStarted: () => {
        clearChallengePolling();
        
        // Find the original tour step
        const stepIndex = driverRef.current?.getActiveIndex() ?? 0;
        const originalStep = tourSteps[stepIndex];
        
        // Save progress so user can resume if they close the tour
        saveTourProgress(tourId, stepIndex);
        
        // Focus panel if specified (brings tab to front)
        if (originalStep?.focusPanel) {
          openPanel(originalStep.focusPanel);
          // Small delay to let the panel render before driver.js positions the highlight
          // Driver.js will automatically reposition after this
        }
        
        // Run onBeforeHighlight callback if provided
        if (originalStep?.onBeforeHighlight) {
          originalStep.onBeforeHighlight();
        }
        
        // Handle challenge if present
        if (originalStep?.challengeId) {
          const challenge = getChallenge(originalStep.challengeId);
          if (challenge) {
            currentChallengeRef.current = challenge;
            
            // Track previous state for change detection
            let lastPanelVisible = challenge.requiredPanel ? isPanelVisible(challenge.requiredPanel) : true;
            let lastComplete = challenge.validate();
            
            // Start polling for challenge completion and panel visibility
            challengeIntervalRef.current = setInterval(() => {
              const panelVisible = challenge.requiredPanel ? isPanelVisible(challenge.requiredPanel) : true;
              const isComplete = challenge.validate();
              
              // Update UI if state changed
              if (panelVisible !== lastPanelVisible || isComplete !== lastComplete) {
                lastPanelVisible = panelVisible;
                lastComplete = isComplete;
                
                // Update the popover content
                const popoverEl = document.querySelector('.driver-popover-description');
                if (popoverEl && originalStep.popover) {
                  const baseDesc = originalStep.popover.description ?? '';
                  popoverEl.innerHTML = `${baseDesc}${buildChallengeHTML(challenge, isComplete)}`;
                }
                
                // Enable/disable next button based on completion
                const nextBtn = document.querySelector('.driver-popover-next-btn') as HTMLButtonElement;
                if (nextBtn) {
                  if (isComplete) {
                    nextBtn.disabled = false;
                    nextBtn.classList.remove('driver-popover-btn--disabled');
                  } else {
                    nextBtn.disabled = true;
                    nextBtn.classList.add('driver-popover-btn--disabled');
                  }
                }
              }
              
              // Stop polling once complete
              if (isComplete) {
                clearChallengePolling();
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
            // Add skip button (unless noSkip is set) and disable next
            const footer = popover.footerButtons;
            const nextBtn = footer.querySelector('.driver-popover-next-btn') as HTMLButtonElement;
            
            if (nextBtn) {
              // Only add skip button if challenge allows skipping
              if (!challenge.noSkip) {
                const skipBtn = document.createElement('button');
                skipBtn.className = 'driver-popover-skip-btn';
                skipBtn.textContent = 'Skip';
                skipBtn.onclick = () => {
                  clearChallengePolling();
                  driverRef.current?.moveNext();
                };
                
                // Insert skip before next
                nextBtn.parentNode?.insertBefore(skipBtn, nextBtn);
              }
              
              // Disable next until challenge complete
              nextBtn.disabled = true;
              nextBtn.classList.add('driver-popover-btn--disabled');
            }
          }
        }
      },

      onDestroyed: () => {
        clearChallengePolling();
        // Check if we completed all steps (progress saved at tourSteps.length - 1 means we were on last step)
        const savedProgress = getTourProgress(tourId);
        if (savedProgress >= tourSteps.length - 1) {
          markTourComplete(tourId);
          clearTourProgress(tourId);
        }
        setIsActive(false);
        setCurrentTour(null);
      },

      onCloseClick: () => {
        // Save progress when user closes mid-tour
        const currentStepIndex = driverRef.current?.getActiveIndex() ?? 0;
        saveTourProgress(tourId, currentStepIndex);
        // Must explicitly destroy - driver.js doesn't do it automatically
        driverRef.current?.destroy();
      },
    });

    return driverRef.current;
  }, [clearChallengePolling, buildChallengeHTML, openPanel]);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      clearChallengePolling();
      if (driverRef.current) {
        driverRef.current.destroy();
      }
    };
  }, [clearChallengePolling]);

  const startTour = useCallback((tourId: TourId, fromBeginning = false) => {
    const tourSteps = tours[tourId];
    if (!tourSteps) {
      console.warn(`Tour "${tourId}" not found`);
      return;
    }

    // If starting from beginning, clear any saved progress
    if (fromBeginning) {
      clearTourProgress(tourId);
    }

    setCurrentTour(tourId);
    setIsActive(true);

    // Small delay to ensure DOM is ready
    setTimeout(() => {
      const driverInstance = getDriverInstance(tourSteps, tourId);
      const driverSteps = convertToDriverSteps(tourSteps);
      driverInstance.setSteps(driverSteps);
      
      // Resume from saved progress (if any)
      const savedStep = fromBeginning ? 0 : getTourProgress(tourId);
      const startStep = Math.min(savedStep, tourSteps.length - 1);
      
      driverInstance.drive(startStep);
    }, 100);
  }, [getDriverInstance, convertToDriverSteps]);

  const endTour = useCallback(() => {
    clearChallengePolling();
    if (driverRef.current) {
      // Save progress (don't mark complete unless they finished)
      if (currentTour) {
        const currentStepIndex = driverRef.current.getActiveIndex() ?? 0;
        saveTourProgress(currentTour, currentStepIndex);
      }
      driverRef.current.destroy();
    }
    setIsActive(false);
    setCurrentTour(null);
  }, [currentTour, clearChallengePolling]);

  const hasCompletedTour = useCallback((tourId: TourId): boolean => {
    return checkTourCompleted(tourId);
  }, []);

  // Check if a tour has been started (has progress or is complete)
  // Useful for preventing auto-start on tours user has already seen
  const hasStartedTour = useCallback((tourId: TourId): boolean => {
    return checkTourCompleted(tourId) || getTourProgress(tourId) > 0;
  }, []);

  const resetAllTours = useCallback(() => {
    resetOnboarding();
  }, []);

  const contextValue: OnboardingContextType = {
    startTour,
    endTour,
    hasCompletedTour,
    hasStartedTour,
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
