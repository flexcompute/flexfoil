import { useEffect, useState } from 'react';
import { DockingLayout } from './components/DockingLayout';
import { ThemeProvider } from './contexts/ThemeContext';
import { OnboardingProvider, useOnboarding } from './onboarding';
import { initWasm } from './lib/wasm';
import { useAirfoilStore } from './stores/airfoilStore';
import { useRunStore } from './stores/runStore';
import { useVisualizationStore } from './stores/visualizationStore';
import { loadFromUrl } from './lib/urlState';
import 'flexlayout-react/style/light.css';
import './App.css';

function AppContent() {
  const [wasmStatus, setWasmStatus] = useState<'loading' | 'ready' | 'error'>('loading');
  const [initialViewport, setInitialViewport] = useState<{ centerX: number; centerY: number; zoom: number } | null>(null);
  const initializeDefaultAirfoil = useAirfoilStore((s) => s.initializeDefaultAirfoil);
  const setControlMode = useAirfoilStore((s) => s.setControlMode);
  const setNPanels = useAirfoilStore((s) => s.setNPanels);
  const setSpacingKnots = useAirfoilStore((s) => s.setSpacingKnots);
  const setDisplayAlpha = useAirfoilStore((s) => s.setDisplayAlpha);
  const setThicknessScale = useAirfoilStore((s) => s.setThicknessScale);
  const setCamberScale = useAirfoilStore((s) => s.setCamberScale);
  const setCurvatureWeight = useAirfoilStore((s) => s.setCurvatureWeight);
  const setSpacingPanelMode = useAirfoilStore((s) => s.setSpacingPanelMode);
  const setSSPInterpolation = useAirfoilStore((s) => s.setSSPInterpolation);
  const setSSPVisualization = useAirfoilStore((s) => s.setSSPVisualization);
  const generateNaca4 = useAirfoilStore((s) => s.generateNaca4);
  
  // Get visualization store setters (individual selectors for stable references)
  const setShowGrid = useVisualizationStore((s) => s.setShowGrid);
  const setShowCurve = useVisualizationStore((s) => s.setShowCurve);
  const setShowPanels = useVisualizationStore((s) => s.setShowPanels);
  const setShowPoints = useVisualizationStore((s) => s.setShowPoints);
  const setShowControls = useVisualizationStore((s) => s.setShowControls);
  const setShowStreamlines = useVisualizationStore((s) => s.setShowStreamlines);
  const setShowSmoke = useVisualizationStore((s) => s.setShowSmoke);
  const setShowPsiContours = useVisualizationStore((s) => s.setShowPsiContours);
  const setShowCp = useVisualizationStore((s) => s.setShowCp);
  const setShowForces = useVisualizationStore((s) => s.setShowForces);
  const setEnableMorphing = useVisualizationStore((s) => s.setEnableMorphing);
  const setMorphDuration = useVisualizationStore((s) => s.setMorphDuration);
  const setStreamlineDensity = useVisualizationStore((s) => s.setStreamlineDensity);
  const setAdaptiveStreamlines = useVisualizationStore((s) => s.setAdaptiveStreamlines);
  const setSmokeDensity = useVisualizationStore((s) => s.setSmokeDensity);
  const setSmokeParticlesPerBlob = useVisualizationStore((s) => s.setSmokeParticlesPerBlob);
  const setFlowSpeed = useVisualizationStore((s) => s.setFlowSpeed);
  const setCpDisplayMode = useVisualizationStore((s) => s.setCpDisplayMode);
  const setCpBarScale = useVisualizationStore((s) => s.setCpBarScale);
  const setForceScale = useVisualizationStore((s) => s.setForceScale);
  
  const initRunDb = useRunStore((s) => s.init);
  const { startTour, hasStartedTour } = useOnboarding();

  // Initialize WASM, run DB, and hydrate from URL
  useEffect(() => {
    initRunDb().catch(err => console.warn('Run DB init failed:', err));
    initWasm()
      .then(() => {
        setWasmStatus('ready');
        
        // Load URL state first
        const urlState = loadFromUrl();
        
        if (urlState) {
          // Apply airfoil settings
          if (urlState.naca && !urlState.custom) {
            const designation = parseInt(urlState.naca, 10);
            const m = Math.floor(designation / 1000) / 100;
            const p = Math.floor((designation % 1000) / 100) / 10;
            const t = (designation % 100) / 100;
            generateNaca4({ m, p, t, nPoints: 100 });
          } else {
            initializeDefaultAirfoil();
          }
          
          if (urlState.nPanels) setNPanels(urlState.nPanels);
          if (urlState.spacing) setSpacingKnots(urlState.spacing);
          if (urlState.mode) setControlMode(urlState.mode);
          if (urlState.alpha !== undefined) setDisplayAlpha(urlState.alpha);
          if (urlState.thicknessScale !== undefined) setThicknessScale(urlState.thicknessScale);
          if (urlState.camberScale !== undefined) setCamberScale(urlState.camberScale);
          if (urlState.curvatureWeight !== undefined) setCurvatureWeight(urlState.curvatureWeight);
          if (urlState.spacingPanelMode) setSpacingPanelMode(urlState.spacingPanelMode);
          if (urlState.sspInterpolation) setSSPInterpolation(urlState.sspInterpolation);
          if (urlState.sspVisualization) setSSPVisualization(urlState.sspVisualization);
          
          // Apply visualization settings
          if (urlState.showGrid !== undefined) setShowGrid(urlState.showGrid);
          if (urlState.showCurve !== undefined) setShowCurve(urlState.showCurve);
          if (urlState.showPanels !== undefined) setShowPanels(urlState.showPanels);
          if (urlState.showPoints !== undefined) setShowPoints(urlState.showPoints);
          if (urlState.showControls !== undefined) setShowControls(urlState.showControls);
          if (urlState.showStreamlines !== undefined) setShowStreamlines(urlState.showStreamlines);
          if (urlState.showSmoke !== undefined) setShowSmoke(urlState.showSmoke);
          if (urlState.showPsiContours !== undefined) setShowPsiContours(urlState.showPsiContours);
          if (urlState.showCp !== undefined) setShowCp(urlState.showCp);
          if (urlState.showForces !== undefined) setShowForces(urlState.showForces);
          
          if (urlState.enableMorphing !== undefined) setEnableMorphing(urlState.enableMorphing);
          if (urlState.morphDuration !== undefined) setMorphDuration(urlState.morphDuration);
          
          if (urlState.streamlineDensity !== undefined) setStreamlineDensity(urlState.streamlineDensity);
          if (urlState.adaptiveStreamlines !== undefined) setAdaptiveStreamlines(urlState.adaptiveStreamlines);
          
          if (urlState.smokeDensity !== undefined) setSmokeDensity(urlState.smokeDensity);
          if (urlState.smokeParticlesPerBlob !== undefined) setSmokeParticlesPerBlob(urlState.smokeParticlesPerBlob);
          
          if (urlState.flowSpeed !== undefined) setFlowSpeed(urlState.flowSpeed);
          
          if (urlState.cpDisplayMode !== undefined) setCpDisplayMode(urlState.cpDisplayMode);
          if (urlState.cpBarScale !== undefined) setCpBarScale(urlState.cpBarScale);
          
          if (urlState.forceScale !== undefined) setForceScale(urlState.forceScale);
          
          // Store viewport for canvas to apply
          if (urlState.viewCenterX !== undefined && urlState.viewCenterY !== undefined && urlState.viewZoom !== undefined) {
            setInitialViewport({
              centerX: urlState.viewCenterX,
              centerY: urlState.viewCenterY,
              zoom: urlState.viewZoom,
            });
          }
        } else {
          // No URL state, initialize defaults
          initializeDefaultAirfoil();
        }
      })
      .catch((err) => {
        console.error('WASM init failed:', err);
        setWasmStatus('error');
      });
  // Note: All store setters are stable functions from zustand, so they don't need to be in deps
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Auto-trigger welcome tour on first visit (only if never started)
  useEffect(() => {
    if (wasmStatus === 'ready' && !hasStartedTour('welcome')) {
      // Delay to let the layout render fully
      const timer = setTimeout(() => {
        startTour('welcome');
      }, 800);
      return () => clearTimeout(timer);
    }
  }, [wasmStatus, hasStartedTour, startTour]);

  return (
    <div className="app-container">
      <DockingLayout wasmStatus={wasmStatus} initialViewport={initialViewport} />
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
