import React from 'react';
import '../../flexfoil-ui/src/App.css';
import {createRoot} from 'react-dom/client';
import {FeedbackWidget} from '../../flexfoil-ui/src/components/FeedbackWidget';
import {CookieConsent} from '../../flexfoil-ui/src/components/CookieConsent';
import {SolvePanel} from '../../flexfoil-ui/src/components/panels/SolvePanel';
import {useRunStore} from '../../flexfoil-ui/src/stores/runStore';
import {useRouteUiStore} from '../../flexfoil-ui/src/stores/routeUiStore';
import {useAirfoilStore} from '../../flexfoil-ui/src/stores/airfoilStore';
const params = new URLSearchParams(location.search);
const mode = params.get('mode');
const consent = params.get('consent') ?? 'granted';
if (consent !== 'missing') localStorage.setItem('ff_cookie_consent', consent);
(window as any).__events = [];
(window as any).__solves = [];
window.gtag = (...args) => (window as any).__events.push(args);
useRunStore.setState({ready: true, lookup: () => null, addRun: async () => {}, addRunBatch: async () => {}, hashPanels: async () => 'offline-fixture'});
useAirfoilStore.setState({solverMode: 'inviscid'});
useRouteUiStore.setState({solveRunMode: mode === 'single_cl' ? 'cl' : 'alpha', solveTargetCl: 0.5,
  solvePolarStart: 0, solvePolarEnd: 1, solvePolarStep: 1,
  sweepPrimary: {param: 'alpha', start: 0, end: 1, step: 1},
  sweepSecondary: mode === 'sweep_2d' ? {param: 'mach', start: 0, end: 0.1, step: 0.1} : null});
createRoot(document.getElementById('root')!).render(<><SolvePanel /><FeedbackWidget /><CookieConsent /><button onClick={() => window.dispatchEvent(new Event('flexfoil:analytics-preferences'))}>Analytics preferences</button></>);
