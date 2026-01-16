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

// Re-export the RustFoil class for advanced usage
export { RustFoil };
