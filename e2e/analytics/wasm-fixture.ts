// An offline solver boundary for UI dispatch tests; these values are not CFD evidence.
export const isWasmReady = () => true;
export function analyzeAirfoilInviscid(_panels: unknown, alpha: number) {
  const result = {success: true, converged: true, cl: alpha * 0.1, cd: 0.01, cm: 0, cp: [], cp_x: [], gamma: [],
    psi_0: 0, iterations: 1, residual: 0, x_tr_upper: 1, x_tr_lower: 1};
  (window as any).__solves.push({alpha, cl: result.cl});
  return result;
}
export const analyzeAirfoil = analyzeAirfoilInviscid;
const unused = () => {throw new Error('Unexpected geometry call in offline analytics fixture');};
export const generateNaca4 = unused, generateNaca4Xfoil = unused,
  repanelWithSpacingAndCurvature = unused, repanelXfoil = unused,
  deflectFlap = unused;
