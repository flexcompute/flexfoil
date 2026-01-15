import { useEffect, useState } from 'react';
import { DockingLayout } from './components/DockingLayout';
import { initWasm } from './lib/wasm';
import 'flexlayout-react/style/dark.css';
import './App.css';

function App() {
  const [wasmStatus, setWasmStatus] = useState<'loading' | 'ready' | 'error'>('loading');

  useEffect(() => {
    initWasm()
      .then(() => setWasmStatus('ready'))
      .catch((err) => {
        console.error('WASM init failed:', err);
        setWasmStatus('error');
      });
  }, []);

  return (
    <div className="app-container">
      <header className="app-header">
        <h1>FlexFoil</h1>
        <div className="wasm-status">
          {wasmStatus === 'loading' && <span className="status-loading">⏳ Loading WASM...</span>}
          {wasmStatus === 'ready' && <span className="status-ready">✓ WASM Ready</span>}
          {wasmStatus === 'error' && <span className="status-error">⚠ WASM Error (using JS fallback)</span>}
        </div>
      </header>
      <main className="app-main">
        <DockingLayout />
      </main>
    </div>
  );
}

export default App;
