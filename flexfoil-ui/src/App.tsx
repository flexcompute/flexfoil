import { useEffect, useState } from 'react';
import { DockingLayout } from './components/DockingLayout';
import { ThemeProvider } from './contexts/ThemeContext';
import { initWasm } from './lib/wasm';
import { hydrateFromUrl, subscribeToUrlSync, useAirfoilStore } from './stores/airfoilStore';
import 'flexlayout-react/style/light.css';
import './App.css';

function AppContent() {
  const [wasmStatus, setWasmStatus] = useState<'loading' | 'ready' | 'error'>('loading');
  const initializeDefaultAirfoil = useAirfoilStore((s) => s.initializeDefaultAirfoil);

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

  return (
    <div className="app-container">
      <DockingLayout wasmStatus={wasmStatus} />
    </div>
  );
}

function App() {
  return (
    <ThemeProvider>
      <AppContent />
    </ThemeProvider>
  );
}

export default App;
