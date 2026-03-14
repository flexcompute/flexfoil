import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import { ModuleRegistry, AllCommunityModule } from 'ag-grid-community';
import {
  AllEnterpriseModule,
  IntegratedChartsModule,
  LicenseManager,
} from 'ag-grid-enterprise';
import { AgChartsEnterpriseModule } from 'ag-charts-enterprise';
import App from './App.tsx'

ModuleRegistry.registerModules([
  AllCommunityModule,
  AllEnterpriseModule,
  IntegratedChartsModule.with(AgChartsEnterpriseModule),
]);

const agLicenseKey = import.meta.env.VITE_AG_GRID_LICENSE_KEY;
if (agLicenseKey) {
  LicenseManager.setLicenseKey(agLicenseKey);
}

createRoot(document.getElementById('root')!).render(
  <StrictMode>
    <App />
  </StrictMode>,
)
