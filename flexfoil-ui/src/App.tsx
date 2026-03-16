import { useCallback, useEffect, useRef, useState } from 'react';
import { useShallow } from 'zustand/react/shallow';
import { DockingLayout } from './components/DockingLayout';
import { ThemeProvider, useTheme } from './contexts/ThemeContext';
import { OnboardingProvider, useOnboarding } from './onboarding';
import { initWasm } from './lib/wasm';
import { useAirfoilStore } from './stores/airfoilStore';
import { DEFAULT_ROUTE_UI_STATE, useRouteUiStore } from './stores/routeUiStore';
import { useRunStore } from './stores/runStore';
import { useVisualizationStore } from './stores/visualizationStore';
import { buildRouteStateSnapshot, parseRouteStateFromLocation, writeRouteState } from './lib/routeState';
import 'flexlayout-react/style/light.css';
import './App.css';

function AppContent() {
  const [wasmStatus, setWasmStatus] = useState<'loading' | 'ready' | 'error'>('loading');
  const initializeDefaultAirfoil = useAirfoilStore((s) => s.initializeDefaultAirfoil);
  const generateNaca4 = useAirfoilStore((s) => s.generateNaca4);
  const hydrateRouteState = useAirfoilStore((s) => s.hydrateRouteState);
  const hydrateVisualizationState = useVisualizationStore((s) => s.hydrateVisualizationState);
  const hydrateUiFromRoute = useRouteUiStore((s) => s.hydrateFromRoute);
  const initRunDb = useRunStore((s) => s.init);
  const { startTour, hasStartedTour } = useOnboarding();
  const { theme, setTheme } = useTheme();
  const basePathRef = useRef('/');
  const isApplyingRouteRef = useRef(false);
  const hasHydratedRouteRef = useRef(false);
  const previousPanelRef = useRef<string | null>(null);
  const wasmStatusRef = useRef(wasmStatus);

  useEffect(() => {
    wasmStatusRef.current = wasmStatus;
  }, [wasmStatus]);

  const airfoilRouteState = useAirfoilStore(
    useShallow((state) => ({
      name: state.name,
      coordinates: state.coordinates,
      panels: state.panels,
      baseCoordinates: state.baseCoordinates,
      camberControlPoints: state.camberControlPoints,
      thicknessControlPoints: state.thicknessControlPoints,
      controlMode: state.controlMode,
      displayAlpha: state.displayAlpha,
      reynolds: state.reynolds,
      mach: state.mach,
      ncrit: state.ncrit,
      maxIterations: state.maxIterations,
      solverMode: state.solverMode,
      thicknessScale: state.thicknessScale,
      camberScale: state.camberScale,
      nPanels: state.nPanels,
      spacingKnots: state.spacingKnots,
      curvatureWeight: state.curvatureWeight,
      spacingPanelMode: state.spacingPanelMode,
      sspInterpolation: state.sspInterpolation,
      sspVisualization: state.sspVisualization,
    })),
  );
  const visualizationRouteState = useVisualizationStore(
    useShallow((state) => ({
      showGrid: state.showGrid,
      showCurve: state.showCurve,
      showPanels: state.showPanels,
      showPoints: state.showPoints,
      showControls: state.showControls,
      showStreamlines: state.showStreamlines,
      showSmoke: state.showSmoke,
      showPsiContours: state.showPsiContours,
      showCp: state.showCp,
      showForces: state.showForces,
      showBoundaryLayer: state.showBoundaryLayer,
      showWake: state.showWake,
      showDisplacementThickness: state.showDisplacementThickness,
      enableMorphing: state.enableMorphing,
      morphDuration: state.morphDuration,
      streamlineDensity: state.streamlineDensity,
      adaptiveStreamlines: state.adaptiveStreamlines,
      smokeDensity: state.smokeDensity,
      smokeParticlesPerBlob: state.smokeParticlesPerBlob,
      smokeWaveSpacing: state.smokeWaveSpacing,
      flowSpeed: state.flowSpeed,
      cpDisplayMode: state.cpDisplayMode,
      cpBarScale: state.cpBarScale,
      forceScale: state.forceScale,
      blThicknessScale: state.blThicknessScale,
      useGPU: state.useGPU,
    })),
  );
  const routeUiState = useRouteUiStore(
    useShallow((state) => ({
      libraryNacaCode: state.libraryNacaCode,
      libraryNPoints: state.libraryNPoints,
      spacingLiveUpdate: state.spacingLiveUpdate,
      solveRunMode: state.solveRunMode,
      solveShowAdvanced: state.solveShowAdvanced,
      solveTargetCl: state.solveTargetCl,
      solvePolarStart: state.solvePolarStart,
      solvePolarEnd: state.solvePolarEnd,
      solvePolarStep: state.solvePolarStep,
      polarXAxis: state.polarXAxis,
      polarYAxis: state.polarYAxis,
      plotDataSource: state.plotDataSource,
      plotChartType: state.plotChartType,
      plotXField: state.plotXField,
      plotYField: state.plotYField,
      plotGroupBy: state.plotGroupBy,
      plotXScale: state.plotXScale,
      plotYScale: state.plotYScale,
      dataExplorerView: state.dataExplorerView,
      dataExplorerSplomKeys: state.dataExplorerSplomKeys,
      dataExplorerColorBy: state.dataExplorerColorBy,
      dataExplorerFilterModel: state.dataExplorerFilterModel,
      activePanel: state.activePanel,
      layoutJson: state.layoutJson,
      viewport: state.viewport,
    })),
  );

  const applyRouteSnapshot = useCallback((pathname?: string, search?: string, hash?: string) => {
    const snapshot = parseRouteStateFromLocation({
      pathname: pathname ?? window.location.pathname,
      search: search ?? window.location.search,
      hash: hash ?? window.location.hash,
    });

    if (!snapshot) {
      initializeDefaultAirfoil();
      return;
    }

    basePathRef.current = snapshot.basePath;

    if (snapshot.theme) {
      setTheme(snapshot.theme);
    }

    const exactGeometry = snapshot.airfoil.exactGeometry;
    if (exactGeometry) {
      hydrateRouteState({
        name: exactGeometry.name,
        coordinates: exactGeometry.coordinates,
        panels: exactGeometry.panels,
        baseCoordinates: exactGeometry.baseCoordinates,
        camberControlPoints: exactGeometry.camberControlPoints,
        thicknessControlPoints: exactGeometry.thicknessControlPoints,
      });
    } else if (snapshot.airfoil.nacaCode) {
      const designation = snapshot.airfoil.nacaCode.padStart(4, '0').slice(0, 4);
      const m = parseInt(designation[0], 10) / 100;
      const p = parseInt(designation[1], 10) / 10;
      const t = parseInt(designation.slice(2), 10) / 100;
      generateNaca4({
        m,
        p,
        t,
        nPoints: snapshot.airfoil.nacaPoints ?? DEFAULT_ROUTE_UI_STATE.libraryNPoints,
      });
    } else {
      initializeDefaultAirfoil();
    }

    hydrateRouteState({
      controlMode: snapshot.airfoil.controlMode,
      displayAlpha: snapshot.airfoil.displayAlpha,
      reynolds: snapshot.airfoil.reynolds,
      mach: snapshot.airfoil.mach,
      ncrit: snapshot.airfoil.ncrit,
      maxIterations: snapshot.airfoil.maxIterations,
      solverMode: snapshot.airfoil.solverMode,
      thicknessScale: snapshot.airfoil.thicknessScale,
      camberScale: snapshot.airfoil.camberScale,
      nPanels: snapshot.airfoil.nPanels,
      spacingKnots: snapshot.airfoil.spacingKnots,
      curvatureWeight: snapshot.airfoil.curvatureWeight,
      spacingPanelMode: snapshot.airfoil.spacingPanelMode,
      sspInterpolation: snapshot.airfoil.sspInterpolation,
      sspVisualization: snapshot.airfoil.sspVisualization,
    });
    hydrateVisualizationState(snapshot.visualization);
    hydrateUiFromRoute({
      ...snapshot.ui,
      activePanel: snapshot.panel,
      layoutJson: snapshot.layoutJson ?? undefined,
      viewport: snapshot.ui.viewport ?? undefined,
    });
  }, [generateNaca4, hydrateRouteState, hydrateUiFromRoute, hydrateVisualizationState, initializeDefaultAirfoil, setTheme]);

  // Initialize WASM, run DB, and hydrate from URL
  useEffect(() => {
    const handlePopState = () => {
      if (wasmStatusRef.current !== 'ready') return;
      isApplyingRouteRef.current = true;
      applyRouteSnapshot();
      window.setTimeout(() => {
        isApplyingRouteRef.current = false;
      }, 0);
    };

    window.addEventListener('popstate', handlePopState);
    initRunDb().catch(err => console.warn('Run DB init failed:', err));
    initWasm()
      .then(() => {
        setWasmStatus('ready');
        isApplyingRouteRef.current = true;
        applyRouteSnapshot();
        window.setTimeout(() => {
          isApplyingRouteRef.current = false;
          hasHydratedRouteRef.current = true;
        }, 0);
      })
      .catch((err) => {
        console.error('WASM init failed:', err);
        setWasmStatus('error');
      });
    return () => window.removeEventListener('popstate', handlePopState);
  }, [applyRouteSnapshot, initRunDb]);

  useEffect(() => {
    if (wasmStatus !== 'ready' || isApplyingRouteRef.current || !hasHydratedRouteRef.current) {
      return;
    }

    const timeoutId = window.setTimeout(() => {
      const snapshot = buildRouteStateSnapshot({
        theme,
        airfoil: airfoilRouteState,
        visualization: visualizationRouteState,
        ui: routeUiState,
      });
      const mode =
        previousPanelRef.current && previousPanelRef.current !== snapshot.panel
          ? 'push'
          : 'replace';
      previousPanelRef.current = snapshot.panel;
      writeRouteState(snapshot, {
        basePath: basePathRef.current,
        mode,
      });
    }, 300);

    return () => window.clearTimeout(timeoutId);
  }, [airfoilRouteState, routeUiState, theme, visualizationRouteState, wasmStatus]);

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
