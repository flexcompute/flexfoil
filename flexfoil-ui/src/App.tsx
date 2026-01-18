import { useEffect, useState } from 'react';
import { DockingLayout } from './components/DockingLayout';
import { ThemeProvider } from './contexts/ThemeContext';
import { initWasm } from './lib/wasm';
import { hydrateFromUrl, subscribeToUrlSync } from './stores/airfoilStore';
import 'flexlayout-react/style/dark.css';
import './App.css';

function AppContent() {
  const [wasmStatus, setWasmStatus] = useState<'loading' | 'ready' | 'error'>('loading');

  // Initialize WASM
  useEffect(() => {
    initWasm()
      .then(() => {
        setWasmStatus('ready');
        // Hydrate from URL after WASM is ready
        hydrateFromUrl();
      })
      .catch((err) => {
        console.error('WASM init failed:', err);
        setWasmStatus('error');
      });
  }, []);

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
