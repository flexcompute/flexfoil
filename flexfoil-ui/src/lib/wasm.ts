/**
 * WASM module initialization and utilities for RustFoil.
 */

import init, {
    generate_naca4,
    generate_naca4_from_string,
    repanel_with_spacing,
    repanel_with_spacing_and_curvature,
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
