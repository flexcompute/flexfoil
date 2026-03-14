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
    analyze_airfoil_faithful,
    compute_streamlines_faithful,
    compute_psi_grid_faithful,
    get_bl_distribution_faithful,
    get_bl_visualization_faithful,
    WasmSmokeSystem,
    greet,
    RustFoil,
    // Morphing interpolation functions
    lerp_points,
    lerp_array,
    lerp_streamlines,
    lerp_morph_state,
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
    cd: number;
    cm: number;
    cp: number[];
    cp_x: number[];
    gamma: number[];  // Now always returned from solver
    psi_0: number;    // Dividing streamline value
    converged: boolean;
    iterations: number;
    residual: number;
    x_tr_upper: number;
    x_tr_lower: number;
    success: boolean;
    error?: string;
}

export interface BLDistribution {
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
    converged: boolean;
    iterations: number;
    residual: number;
    success: boolean;
    error?: string;
}

export interface BLStationData {
    x: number[];
    y: number[];
    delta_star: number[];
    theta: number[];
    cf: number[];
    ue: number[];
    hk: number[];
    ampl: number[];
    is_laminar: boolean[];
}

export interface WakeData {
    x: number[];
    y: number[];
    delta_star: number[];
    theta: number[];
}

export interface BLVisualizationData {
    upper: BLStationData;
    lower: BLStationData;
    wake: WakeData;
    wake_upper_fraction: number;
    wake_geometry_x: number[];
    wake_geometry_y: number[];
    x_tr_upper: number;
    x_tr_lower: number;
    converged: boolean;
    iterations: number;
    residual: number;
    success: boolean;
    error?: string;
}

export function getBLVisualizationData(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number = 1e6,
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100
): BLVisualizationData {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const coordsFlat = pointsToFlat(coordinates);
    return get_bl_visualization_faithful(
        coordsFlat,
        alphaDeg,
        reynolds,
        mach,
        ncrit,
        maxIterations
    ) as BLVisualizationData;
}

export function analyzeAirfoil(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number = 1e6,
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100
): AnalysisResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return analyze_airfoil_faithful(
        coordsFlat,
        alphaDeg,
        reynolds,
        mach,
        ncrit,
        maxIterations
    ) as AnalysisResult;
}

export function getBLDistribution(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number = 1e6,
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100
): BLDistribution {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const coordsFlat = pointsToFlat(coordinates);
    return get_bl_distribution_faithful(
        coordsFlat,
        alphaDeg,
        reynolds,
        mach,
        ncrit,
        maxIterations
    ) as BLDistribution;
}

/**
 * Compute gamma (circulation) values for each panel node.
 * 
 * This calls the actual panel method solver to get the correct gamma values
 * needed for velocity field computation (GPU visualization).
 * Returns an array of gamma values at each panel node.
 */
export interface GammaResult {
    gamma: Float64Array | null;
    psi_0: number;
    success: boolean;
    error?: string;
}

export function computeGamma(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number = 1e6,
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100
): GammaResult {
    if (!initialized) {
        return { gamma: null, psi_0: 0, success: false, error: 'WASM not initialized' };
    }
    
    try {
        // Call the actual solver to get real gamma values
        const result = analyzeAirfoil(coordinates, alphaDeg, reynolds, mach, ncrit, maxIterations);
        
        if (!result.success || !result.gamma || result.gamma.length === 0) {
            return { 
                gamma: null, 
                psi_0: 0, 
                success: false, 
                error: result.error || 'Solver returned no gamma values' 
            };
        }
        
        // Convert to Float64Array
        const gamma = new Float64Array(result.gamma);
        
        return { 
            gamma, 
            psi_0: result.psi_0 || 0,
            success: true 
        };
    } catch (e) {
        return { gamma: null, psi_0: 0, success: false, error: String(e) };
    }
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
    reynolds: number = 1e6,
    seedCount: number = 25,
    bounds: [number, number, number, number] = [-0.5, 2.0, -0.5, 0.5],
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100
): StreamlineResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return compute_streamlines_faithful(
        coordsFlat,
        alphaDeg,
        reynolds,
        mach,
        ncrit,
        maxIterations,
        seedCount,
        new Float64Array(bounds)
    ) as StreamlineResult;
}

/**
 * Stream function grid result from WASM.
 */
export interface PsiGridResult {
    /** Stream function values in row-major order (ny rows, nx columns) */
    grid: number[];
    /** Internal stream function value (contour at this value = dividing streamline) */
    psi_0: number;
    /** Grid width */
    nx: number;
    /** Grid height */
    ny: number;
    /** Minimum psi value (excluding interior NaN) */
    psi_min: number;
    /** Maximum psi value (excluding interior NaN) */
    psi_max: number;
    success: boolean;
    error?: string;
}

/**
 * Compute stream function values on a grid.
 * 
 * The stream function ψ is constant along streamlines. The contour ψ = psi_0
 * represents the dividing streamline that passes through the stagnation point.
 * 
 * @param coordinates - Airfoil coordinates
 * @param alphaDeg - Angle of attack in degrees
 * @param bounds - [xMin, xMax, yMin, yMax] for grid domain
 * @param resolution - [nx, ny] grid resolution
 */
export function computePsiGrid(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number = 1e6,
    bounds: [number, number, number, number] = [-1.0, 2.0, -1.0, 1.0],
    resolution: [number, number] = [100, 80],
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100
): PsiGridResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    return compute_psi_grid_faithful(
        coordsFlat, 
        alphaDeg, 
        reynolds,
        mach,
        ncrit,
        maxIterations,
        new Float64Array(bounds), 
        new Uint32Array(resolution)
    ) as PsiGridResult;
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
// Morphing Interpolation (WASM-accelerated)
// ============================================================================

/**
 * Interpolate between two point arrays using WASM.
 * Much faster than JavaScript for large arrays.
 */
export function wasmLerpPoints(
    from: { x: number; y: number }[],
    to: { x: number; y: number }[],
    t: number
): { x: number; y: number }[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const fromFlat = pointsToFlat(from);
    const toFlat = pointsToFlat(to);
    const result = lerp_points(fromFlat, toFlat, t);
    return flatToPoints(result);
}

/**
 * Interpolate between two number arrays using WASM.
 */
export function wasmLerpArray(
    from: number[],
    to: number[],
    t: number
): number[] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const result = lerp_array(new Float64Array(from), new Float64Array(to), t);
    return Array.from(result);
}

/**
 * Encode streamlines to flat array format for WASM.
 * Format: [n_lines, n_pts_0, x0, y0, x1, y1, ..., n_pts_1, x0, y0, ...]
 */
export function encodeStreamlines(streamlines: [number, number][][]): Float64Array {
    // Calculate total size
    const totalPts = streamlines.reduce((sum, line) => sum + line.length, 0);
    const totalSize = 1 + streamlines.length + totalPts * 2;
    
    const result = new Float64Array(totalSize);
    result[0] = streamlines.length;
    
    let idx = 1;
    for (const line of streamlines) {
        result[idx++] = line.length;
        for (const [x, y] of line) {
            result[idx++] = x;
            result[idx++] = y;
        }
    }
    
    return result;
}

/**
 * Decode streamlines from flat array format.
 */
export function decodeStreamlines(data: Float64Array | number[]): [number, number][][] {
    if (data.length === 0) {
        return [];
    }
    
    const nLines = data[0];
    const result: [number, number][][] = [];
    let idx = 1;
    
    for (let i = 0; i < nLines; i++) {
        if (idx >= data.length) break;
        
        const nPts = data[idx++];
        const line: [number, number][] = [];
        
        for (let j = 0; j < nPts; j++) {
            if (idx + 1 >= data.length) break;
            line.push([data[idx], data[idx + 1]]);
            idx += 2;
        }
        result.push(line);
    }
    
    return result;
}

/**
 * Interpolate between two streamline arrays using WASM.
 */
export function wasmLerpStreamlines(
    from: [number, number][][],
    to: [number, number][][],
    t: number
): [number, number][][] {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const fromEncoded = encodeStreamlines(from);
    const toEncoded = encodeStreamlines(to);
    const result = lerp_streamlines(fromEncoded, toEncoded, t);
    return decodeStreamlines(result);
}

/**
 * Result of WASM morph interpolation.
 */
export interface MorphResult {
    coordinates: number[];
    panels: number[];
    cp: number[];
    cp_x: number[];
}

/**
 * Batch interpolation for morphing - interpolates coordinates, panels, cp, cpX
 * in a single WASM call for maximum performance.
 */
export function wasmLerpMorphState(
    fromCoords: { x: number; y: number }[],
    toCoords: { x: number; y: number }[],
    fromPanels: { x: number; y: number }[],
    toPanels: { x: number; y: number }[],
    fromCp: number[],
    toCp: number[],
    fromCpX: number[],
    toCpX: number[],
    t: number
): {
    coordinates: { x: number; y: number }[];
    panels: { x: number; y: number }[];
    cp: number[];
    cpX: number[];
} {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const result = lerp_morph_state(
        pointsToFlat(fromCoords),
        pointsToFlat(toCoords),
        pointsToFlat(fromPanels),
        pointsToFlat(toPanels),
        new Float64Array(fromCp),
        new Float64Array(toCp),
        new Float64Array(fromCpX),
        new Float64Array(toCpX),
        t
    ) as MorphResult;
    
    return {
        coordinates: flatToPoints(result.coordinates),
        panels: flatToPoints(result.panels),
        cp: Array.from(result.cp),
        cpX: Array.from(result.cp_x),
    };
}

// Re-export the RustFoil class and WasmSmokeSystem for advanced usage
export { RustFoil, WasmSmokeSystem };
