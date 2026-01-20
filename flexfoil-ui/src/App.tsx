import { useEffect, useState } from 'react';
import { DockingLayout } from './components/DockingLayout';
import { ThemeProvider } from './contexts/ThemeContext';
import { OnboardingProvider, useOnboarding } from './onboarding';
import { initWasm } from './lib/wasm';
import { hydrateFromUrl, subscribeToUrlSync, useAirfoilStore } from './stores/airfoilStore';
import 'flexlayout-react/style/light.css';
import './App.css';

function AppContent() {
  const [wasmStatus, setWasmStatus] = useState<'loading' | 'ready' | 'error'>('loading');
  const initializeDefaultAirfoil = useAirfoilStore((s) => s.initializeDefaultAirfoil);
  const { startTour, hasCompletedTour } = useOnboarding();

  // Initialize WASM
  useEffect(() => {
    initWasm()
      .then(() => {
        setWasmStatus('ready');
        // Initialize default NACA 0012 with proper geometry
        initializeDefaultAirfoil();
        // Hydrate from URL after WASM is ready (may override default)
        hydrateFromUrl();
      })
      .catch((err) => {
        console.error('WASM init failed:', err);
        setWasmStatus('error');
      });
  }, [initializeDefaultAirfoil]);

  // Subscribe to URL sync
  useEffect(() => {
    const unsubscribe = subscribeToUrlSync();
    return () => unsubscribe();
  }, []);

  // Auto-trigger welcome tour on first visit
  useEffect(() => {
    if (wasmStatus === 'ready' && !hasCompletedTour('welcome')) {
      // Delay to let the layout render fully
      const timer = setTimeout(() => {
        startTour('welcome');
      }, 800);
      return () => clearTimeout(timer);
    }
  }, [wasmStatus, hasCompletedTour, startTour]);

  return (
    <div className="app-container">
      <DockingLayout wasmStatus={wasmStatus} />
    </div>
  );
}

function App() {
  return (
    <ThemeProvider>
      <OnboardingProvider>
        <AppContent />
      </OnboardingProvider>
    </ThemeProvider>
  );
}

export default App;
