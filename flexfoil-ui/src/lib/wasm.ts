/**
 * WASM module initialization and utilities for RustFoil.
 */

import init, {
    generate_naca4,
    generate_naca4_from_string,
    generate_naca4_xfoil,
    repanel_with_spacing_and_curvature,
    repanel_xfoil,
    repanel_xfoil_with_params,
    compute_curvature_spacing,
    analyze_airfoil,
    analyze_viscous,
    analyze_viscous_extended,
    get_bl_distribution,
    get_turbulent_model_info,
    run_re_sweep,
    run_viscous_polar,
    compute_streamlines,
    WasmSmokeSystem,
    greet,
    RustFoil,
} from 'rustfoil-wasm';

let initialized = false;
let initPromise: Promise<void> | null = null;

/**
 * Initialize the WASM module. Safe to call multiple times.
 */
export async function initWasm(): Promise<void> {
    if (initialized) return;
    if (initPromise) return initPromise;
    
    initPromise = (async () => {
        await init();
        initialized = true;
        console.log('RustFoil WASM initialized:', greet());
    })();
    
    return initPromise;
}

/**
 * Check if WASM is initialized.
 */
export function isWasmReady(): boolean {
    return initialized;
}

/**
 * Convert flat coordinate array to Point array.
 */
export function flatToPoints(flat: Float64Array | number[]): { x: number; y: number }[] {
    const points: { x: number; y: number }[] = [];
    for (let i = 0; i < flat.length; i += 2) {
        points.push({ x: flat[i], y: flat[i + 1] });
    }
    return points;
}

/**
 * Convert Point array to flat coordinate array.
 */
export function pointsToFlat(points: { x: number; y: number }[]): Float64Array {
    const flat = new Float64Array(points.length * 2);
    for (let i = 0; i < points.length; i++) {
        flat[i * 2] = points[i].x;
        flat[i * 2 + 1] = points[i].y;
    }
    return flat;
}

/**
 * Generate NACA 4-series airfoil coordinates.
 * 
 * @param m - Maximum camber (0-9, first digit)
 * @param p - Position of max camber (0-9, second digit)
 * @param t - Thickness (00-99, last two digits)
 * @param nPoints - Number of points per surface
 * @returns Array of {x, y} points
 */
export function generateNaca4(
    m: number, 
    p: number, 
    t: number, 
    nPoints: number
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    // Convert NACA digits to fractions
    const mFrac = m / 100;  // First digit / 100
    const pFrac = p / 10;   // Second digit / 10
    const tFrac = t / 100;  // Last two digits / 100
    
    const flat = generate_naca4(mFrac, pFrac, tFrac, nPoints);
    return flatToPoints(flat);
}

/**
 * Generate NACA 4-series from string designation (e.g., "2412", "0012").
 */
export function generateNaca4FromString(
    digits: string, 
    nPoints: number
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const flat = generate_naca4_from_string(digits, nPoints);
    if (flat.length === 0) {
        throw new Error(`Invalid NACA designation: ${digits}`);
    }
    return flatToPoints(flat);
}

/**
 * Generate NACA 4-series using XFOIL's exact algorithm.
 * 
 * This matches XFOIL's naca.f SUBROUTINE NACA4 exactly:
 * - TE bunching parameter AN = 1.5
 * - Blunt trailing edge (not closed)
 * - Perfectly symmetric output for symmetric airfoils
 * - Output order: upper TE → LE → lower TE
 * 
 * @param designation - 4-digit NACA designation as number (e.g., 12 for NACA 0012, 2412 for NACA 2412)
 * @param nPointsPerSide - Number of points per side (default: 123 = XFOIL's IQX/3)
 * @returns Array of {x, y} points. Total points = 2*n - 1.
 */
export function generateNaca4Xfoil(
    designation: number,
    nPointsPerSide?: number
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const flat = generate_naca4_xfoil(designation, nPointsPerSide);
    return flatToPoints(flat);
}

/**
 * Repanel airfoil with custom SSP spacing.
 * 
 * @param coordinates - Current airfoil coordinates
 * @param spacingKnots - Array of {s, f} knots for spacing function
 * @param nPanels - Desired number of panels
 */
export function repanelWithSpacing(
    coordinates: { x: number; y: number }[],
    spacingKnots: { s: number; f: number }[],
    nPanels: number
): { x: number; y: number }[] {
    return repanelWithSpacingAndCurvature(coordinates, spacingKnots, nPanels, 0);
}

/**
 * Repanel airfoil with blended SSP and curvature-based spacing.
 * 
 * @param coordinates - Current airfoil coordinates
 * @param spacingKnots - Array of {s, f} knots for spacing function
 * @param nPanels - Desired number of panels
 * @param curvatureWeight - Blend factor: 0.0 = pure SSP, 1.0 = pure curvature-based
 */
export function repanelWithSpacingAndCurvature(
    coordinates: { x: number; y: number }[],
    spacingKnots: { s: number; f: number }[],
    nPanels: number,
    curvatureWeight: number
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    const knotsFlat = new Float64Array(spacingKnots.length * 2);
    for (let i = 0; i < spacingKnots.length; i++) {
        knotsFlat[i * 2] = spacingKnots[i].s;
        knotsFlat[i * 2 + 1] = spacingKnots[i].f;
    }
    
    const result = repanel_with_spacing_and_curvature(coordsFlat, knotsFlat, nPanels, curvatureWeight);
    return flatToPoints(result);
}

/**
 * Compute curvature-based spacing values along the airfoil.
 * 
 * @param coordinates - Airfoil coordinates
 * @param nSamples - Number of sample points
 * @returns Array of {s, f} values representing curvature-based spacing
 */
export function computeCurvatureSpacing(
    coordinates: { x: number; y: number }[],
    nSamples: number = 100
): { s: number; f: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    const result = compute_curvature_spacing(coordsFlat, nSamples);
    
    // Convert flat array to {s, f} pairs
    const knots: { s: number; f: number }[] = [];
    for (let i = 0; i < result.length; i += 2) {
        knots.push({ s: result[i], f: result[i + 1] });
    }
    return knots;
}

/**
 * Repanel airfoil using XFOIL's curvature-based algorithm.
 * 
 * This implements Mark Drela's PANGEN algorithm from XFOIL, which distributes
 * panel nodes based on local curvature. Panels are bunched in regions of high
 * curvature (leading edge) and at the trailing edge.
 * 
 * @param coordinates - Current airfoil coordinates
 * @param nPanels - Desired number of panels (XFOIL default: 160)
 * @returns Repaneled coordinates
 */
export function repanelXfoil(
    coordinates: { x: number; y: number }[],
    nPanels: number = 160
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    const result = repanel_xfoil(coordsFlat, nPanels);
    return flatToPoints(result);
}

/**
 * XFOIL paneling parameters.
 */
export interface XfoilPanelingParams {
    /** Curvature attraction (0=uniform, 1=XFOIL default bunching) */
    curvParam?: number;
    /** TE/LE panel density ratio (XFOIL default: 0.15) */
    teLeRatio?: number;
    /** TE panel length ratio (XFOIL default: 0.667) */
    teSpacingRatio?: number;
}

/**
 * Repanel airfoil using XFOIL's algorithm with custom parameters.
 * 
 * @param coordinates - Current airfoil coordinates
 * @param nPanels - Desired number of panels
 * @param params - XFOIL paneling parameters
 * @returns Repaneled coordinates
 */
export function repanelXfoilWithParams(
    coordinates: { x: number; y: number }[],
    nPanels: number,
    params: XfoilPanelingParams = {}
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    const result = repanel_xfoil_with_params(
        coordsFlat,
        nPanels,
        params.curvParam ?? 1.0,
        params.teLeRatio ?? 0.15,
        params.teSpacingRatio ?? 0.667
    );
    return flatToPoints(result);
}

/**
 * Analyze airfoil at given angle of attack.
 */
export interface AnalysisResult {
    cl: number;
    cm: number;
    cp: number[];
    cp_x: number[];
    success: boolean;
    error?: string;
}

export function analyzeAirfoil(
    coordinates: { x: number; y: number }[],
    alphaDeg: number
): AnalysisResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return analyze_airfoil(coordsFlat, alphaDeg) as AnalysisResult;
}

/**
 * Streamline result from WASM.
 */
export interface StreamlineResult {
    streamlines: [number, number][][];
    success: boolean;
    error?: string;
}

/**
 * Compute streamlines for flow visualization.
 * 
 * @param coordinates - Airfoil coordinates
 * @param alphaDeg - Angle of attack in degrees
 * @param seedCount - Number of streamlines
 * @param bounds - [xMin, xMax, yMin, yMax] for domain
 */
export function computeStreamlines(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    seedCount: number = 25,
    bounds: [number, number, number, number] = [-0.5, 2.0, -0.5, 0.5]
): StreamlineResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return compute_streamlines(coordsFlat, alphaDeg, seedCount, new Float64Array(bounds)) as StreamlineResult;
}

/**
 * Create a smoke visualization system.
 * 
 * @param spawnY - Y coordinates for spawn points
 * @param spawnX - X coordinate for spawn points
 * @param particlesPerBlob - Particles per blob (default 20)
 */
export function createSmokeSystem(
    spawnY: number[],
    spawnX: number = -0.5,
    particlesPerBlob: number = 20
): WasmSmokeSystem {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    return new WasmSmokeSystem(new Float64Array(spawnY), spawnX, particlesPerBlob);
}

// ============================================================================
// Viscous Analysis
// ============================================================================

/**
 * Turbulent boundary layer model selection.
 * 
 * - Head (0): Head's entrainment method (1958) - fast, good for attached flows
 * - XfoilCtau (1): XFOIL's Cτ lag-dissipation method (Drela 1987) - accurate, good for separation
 * - GreenLag (2): Green's lag-entrainment method (1977) - best history effects
 */
export type TurbulentModel = 0 | 1 | 2;

export const TurbulentModelNames: Record<TurbulentModel, string> = {
    0: 'Head',
    1: 'XFOIL Cτ',
    2: 'Green Lag',
};

export const TurbulentModelDescriptions: Record<TurbulentModel, string> = {
    0: 'Fast, good for attached flows',
    1: 'Accurate, good for separation',
    2: 'Best for history effects',
};

/**
 * Information about available turbulent models.
 */
export interface TurbulentModelInfo {
    models: {
        id: number;
        name: string;
        fullName: string;
        description: string;
        equations: string[];
        recommended: string;
    }[];
    default: string;
}

/**
 * Result of a viscous analysis.
 */
export interface ViscousAnalysisResult {
    cl: number;
    cd: number;
    cd_friction: number;
    cd_pressure: number;
    cm: number;
    cp: number[];
    cp_x: number[];
    x_tr_upper: number;
    x_tr_lower: number;
    converged: boolean;
    iterations: number;
    reynolds: number;
    alpha: number;
    success: boolean;
    error?: string;
}

/**
 * Boundary layer distribution data for plotting.
 */
export interface BLDistribution {
    s_upper: number[];
    s_lower: number[];
    x_upper: number[];
    x_lower: number[];
    theta_upper: number[];
    theta_lower: number[];
    delta_star_upper: number[];
    delta_star_lower: number[];
    h_upper: number[];
    h_lower: number[];
    cf_upper: number[];
    cf_lower: number[];
    x_tr_upper: number;
    x_tr_lower: number;
}

/**
 * Point in a Reynolds number sweep.
 */
export interface ReSweepPoint {
    reynolds: number;
    cl: number;
    cd: number;
    cm: number;
    x_tr_upper: number;
    x_tr_lower: number;
    converged: boolean;
}

/**
 * Viscous polar point.
 */
export interface ViscousPolarPoint {
    alpha: number;
    cl: number;
    cd: number;
    cm: number;
    x_tr_upper: number;
    x_tr_lower: number;
    converged: boolean;
}

/**
 * Analyze airfoil with viscous boundary layer.
 * 
 * @param coordinates - Airfoil coordinates
 * @param alphaDeg - Angle of attack in degrees
 * @param reynolds - Chord Reynolds number
 */
export function analyzeViscous(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number
): ViscousAnalysisResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return analyze_viscous(coordsFlat, alphaDeg, reynolds) as ViscousAnalysisResult;
}

/**
 * Analyze airfoil with viscous boundary layer and extended options.
 * 
 * @param coordinates - Airfoil coordinates
 * @param alphaDeg - Angle of attack in degrees
 * @param reynolds - Chord Reynolds number
 * @param model - Turbulent BL model (0=Head, 1=XfoilCtau, 2=GreenLag)
 * @param nCrit - Critical N-factor for transition (typically 9.0, range 5-12)
 */
export function analyzeViscousExtended(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number,
    model: TurbulentModel,
    nCrit: number
): ViscousAnalysisResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return analyze_viscous_extended(coordsFlat, alphaDeg, reynolds, model, nCrit) as ViscousAnalysisResult;
}

/**
 * Get information about available turbulent models.
 */
export function getTurbulentModelInfo(): TurbulentModelInfo {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    return get_turbulent_model_info() as TurbulentModelInfo;
}

/**
 * Get boundary layer distribution data for plotting.
 * 
 * @param coordinates - Airfoil coordinates
 * @param alphaDeg - Angle of attack in degrees
 * @param reynolds - Chord Reynolds number
 */
export function getBLDistribution(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number
): BLDistribution | null {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return get_bl_distribution(coordsFlat, alphaDeg, reynolds) as BLDistribution | null;
}

/**
 * Run a Reynolds number sweep at fixed alpha.
 * 
 * @param coordinates - Airfoil coordinates
 * @param alphaDeg - Angle of attack in degrees
 * @param reStart - Starting Reynolds number
 * @param reEnd - Ending Reynolds number
 * @param reCount - Number of Reynolds numbers to evaluate
 */
export function runReSweep(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reStart: number,
    reEnd: number,
    reCount: number
): ReSweepPoint[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return run_re_sweep(coordsFlat, alphaDeg, reStart, reEnd, reCount) as ReSweepPoint[];
}

/**
 * Run a viscous alpha polar sweep.
 * 
 * @param coordinates - Airfoil coordinates
 * @param reynolds - Chord Reynolds number
 * @param alphaStart - Starting alpha (degrees)
 * @param alphaEnd - Ending alpha (degrees)
 * @param alphaStep - Alpha step (degrees)
 */
export function runViscousPolar(
    coordinates: { x: number; y: number }[],
    reynolds: number,
    alphaStart: number,
    alphaEnd: number,
    alphaStep: number
): ViscousPolarPoint[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return run_viscous_polar(coordsFlat, reynolds, alphaStart, alphaEnd, alphaStep) as ViscousPolarPoint[];
}

// Re-export the RustFoil class and WasmSmokeSystem for advanced usage
export { RustFoil, WasmSmokeSystem };
