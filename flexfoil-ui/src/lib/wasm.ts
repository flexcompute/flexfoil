/**
 * WASM module initialization and utilities for RustFoil.
 */

import init, {
    generate_naca4,
    generate_naca4_xfoil,
    repanel_with_spacing_and_curvature,
    repanel_xfoil,
    compute_curvature_spacing,
    analyze_airfoil as analyze_airfoil_inviscid_wasm,
    analyze_airfoil_faithful,
    compute_streamlines as compute_streamlines_inviscid_wasm,
    compute_streamlines_faithful,
    compute_dividing_streamline as compute_dividing_streamline_inviscid_wasm,
    compute_dividing_streamline_faithful,
    compute_psi_grid as compute_psi_grid_inviscid_wasm,
    compute_psi_grid_faithful,
    get_bl_distribution_faithful,
    get_bl_visualization_faithful,
    WasmSmokeSystem,
    greet,
    RustFoil,
    lerp_streamlines,
    lerp_morph_state,
    inverse_design_qdes,
    gdes_rotate,
    gdes_scale,
    gdes_translate,
    gdes_set_te_gap,
    gdes_flap,
    gdes_set_le_radius,
    gdes_scale_thickness,
    gdes_scale_camber,
    full_inverse_design_mdes,
    analyze_multi_element as analyze_multi_element_wasm,
    analyze_multi_element_viscous as analyze_multi_element_viscous_wasm,
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

/**
 * Inviscid-only panel method (no boundary layer coupling).
 * Returns CL, CM, Cp — CD is always 0.
 */
export function analyzeAirfoilInviscid(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
): AnalysisResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const coordsFlat = pointsToFlat(coordinates);
    return analyze_airfoil_inviscid_wasm(coordsFlat, alphaDeg) as AnalysisResult;
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

export interface DividingStreamlineResult {
    streamline: [number, number][];
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
    maxIterations: number = 100,
    solverMode: 'inviscid' | 'viscous' = 'viscous'
): StreamlineResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    if (solverMode === 'inviscid') {
        return compute_streamlines_inviscid_wasm(
            coordsFlat,
            alphaDeg,
            seedCount,
            new Float64Array(bounds)
        ) as StreamlineResult;
    }
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

export function computeDividingStreamline(
    coordinates: { x: number; y: number }[],
    alphaDeg: number,
    reynolds: number = 1e6,
    bounds: [number, number, number, number] = [-0.5, 2.0, -0.5, 0.5],
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100,
    solverMode: 'inviscid' | 'viscous' = 'viscous'
): DividingStreamlineResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const coordsFlat = pointsToFlat(coordinates);
    if (solverMode === 'inviscid') {
        return compute_dividing_streamline_inviscid_wasm(
            coordsFlat,
            alphaDeg,
            new Float64Array(bounds)
        ) as DividingStreamlineResult;
    }

    return compute_dividing_streamline_faithful(
        coordsFlat,
        alphaDeg,
        reynolds,
        mach,
        ncrit,
        maxIterations,
        new Float64Array(bounds)
    ) as DividingStreamlineResult;
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
    maxIterations: number = 100,
    solverMode: 'inviscid' | 'viscous' = 'viscous'
): PsiGridResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }
    
    const coordsFlat = pointsToFlat(coordinates);
    if (solverMode === 'inviscid') {
        return compute_psi_grid_inviscid_wasm(
            coordsFlat,
            alphaDeg,
            new Float64Array(bounds),
            new Uint32Array(resolution)
        ) as PsiGridResult;
    }
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

/**
 * Encode streamlines to flat array format for WASM.
 * Format: [n_lines, n_pts_0, x0, y0, x1, y1, ..., n_pts_1, x0, y0, ...]
 */
function encodeStreamlines(streamlines: [number, number][][]): Float64Array {
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
function decodeStreamlines(data: Float64Array | number[]): [number, number][][] {
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

// ============================================================================
// Inverse Design (QDES)
// ============================================================================

export interface InverseDesignIterationSnapshot {
    iteration: number;
    rms_error: number;
    max_error: number;
    geometry_delta_norm: number;
    cl: number;
    cd: number;
}

export interface InverseDesignResult {
    input_coords: number[];
    output_coords: number[];
    paneled_coords: number[];
    cl: number;
    cd: number;
    cm: number;
    converged: boolean;
    iterations: number;
    rms_error: number;
    max_error: number;
    target_kind: 'cp' | 'velocity';
    target_upper_x: number[];
    target_upper_values: number[];
    target_lower_x: number[];
    target_lower_values: number[];
    achieved_upper_x: number[];
    achieved_upper_values: number[];
    achieved_lower_x: number[];
    achieved_lower_values: number[];
    history: InverseDesignIterationSnapshot[];
    success: boolean;
    error?: string;
}

export interface InverseDesignTarget {
    x: number[];
    values: number[];
}

export interface InverseDesignOptions {
    alphaDeg: number;
    reynolds?: number;
    mach?: number;
    ncrit?: number;
    targetKind: 'cp' | 'velocity';
    upper?: InverseDesignTarget;
    lower?: InverseDesignTarget;
    maxDesignIterations?: number;
    damping?: number;
}

/**
 * Run inverse design (QDES) to modify an airfoil to match a target Cp or velocity distribution.
 *
 * Provide target distributions on the upper and/or lower surface.
 * The solver iterates, perturbing the geometry to match the target.
 */
export function runInverseDesign(
    coordinates: { x: number; y: number }[],
    options: InverseDesignOptions,
): InverseDesignResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const coordsFlat = pointsToFlat(coordinates);
    const upperX = options.upper ? new Float64Array(options.upper.x) : new Float64Array(0);
    const upperValues = options.upper ? new Float64Array(options.upper.values) : new Float64Array(0);
    const lowerX = options.lower ? new Float64Array(options.lower.x) : new Float64Array(0);
    const lowerValues = options.lower ? new Float64Array(options.lower.values) : new Float64Array(0);

    return inverse_design_qdes(
        coordsFlat,
        options.alphaDeg,
        options.reynolds ?? 1e6,
        options.mach ?? 0,
        options.ncrit ?? 9,
        options.targetKind,
        upperX,
        upperValues,
        lowerX,
        lowerValues,
        options.maxDesignIterations ?? 6,
        options.damping ?? 0.6,
    ) as InverseDesignResult;
}

// ============================================================================
// Geometry Design (GDES)
// ============================================================================

export interface GeometryResult {
    coords: number[];
    success: boolean;
    error?: string;
}

function gdesResultToPoints(result: GeometryResult): { x: number; y: number }[] {
    if (!result.success || result.coords.length < 2) return [];
    return flatToPoints(result.coords);
}

/** Rotate airfoil by angle (degrees) about center (cx, cy). */
export function rotateAirfoil(
    coordinates: { x: number; y: number }[],
    angleDeg: number,
    cx: number = 0.25,
    cy: number = 0,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_rotate(pointsToFlat(coordinates), angleDeg, cx, cy) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Scale airfoil about center (cx, cy) by factors (sx, sy). */
export function scaleAirfoil(
    coordinates: { x: number; y: number }[],
    sx: number,
    sy: number,
    cx: number = 0,
    cy: number = 0,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_scale(pointsToFlat(coordinates), sx, sy, cx, cy) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Translate airfoil by (dx, dy). */
export function translateAirfoil(
    coordinates: { x: number; y: number }[],
    dx: number,
    dy: number,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_translate(pointsToFlat(coordinates), dx, dy) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Set trailing-edge gap (fraction of chord). blendFraction controls how far forward the change extends. */
export function setTeGap(
    coordinates: { x: number; y: number }[],
    gap: number,
    blendFraction: number = 0.8,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_set_te_gap(pointsToFlat(coordinates), gap, blendFraction) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Apply a trailing-edge flap deflection.
 *  hingeXFrac: hinge x as fraction of chord.
 *  hingeYFrac: hinge y as fraction of local thickness (0=lower, 0.5=mid, 1=upper).
 */
export function deflectFlap(
    coordinates: { x: number; y: number }[],
    hingeXFrac: number,
    hingeYFrac: number,
    deflectionDeg: number,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_flap(pointsToFlat(coordinates), hingeXFrac, hingeYFrac, deflectionDeg) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Scale leading-edge radius by factor (1.0 = unchanged). */
export function setLeRadius(
    coordinates: { x: number; y: number }[],
    factor: number,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_set_le_radius(pointsToFlat(coordinates), factor) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Scale thickness while preserving camber (1.0 = unchanged). */
export function scaleThickness(
    coordinates: { x: number; y: number }[],
    factor: number,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_scale_thickness(pointsToFlat(coordinates), factor) as GeometryResult;
    return gdesResultToPoints(result);
}

/** Scale camber while preserving thickness (1.0 = unchanged, 0.0 = symmetric). */
export function scaleCamber(
    coordinates: { x: number; y: number }[],
    factor: number,
): { x: number; y: number }[] {
    if (!initialized) throw new Error('WASM not initialized.');
    const result = gdes_scale_camber(pointsToFlat(coordinates), factor) as GeometryResult;
    return gdesResultToPoints(result);
}

// ============================================================================
// Full-Inverse Design (MDES) via Circle-Plane Conformal Mapping
// ============================================================================

export interface MdesDesignResult {
    x: number[];
    y: number[];
    qspec_x: number[];
    qspec_values: number[];
    cl: number;
    cm: number;
    success: boolean;
    error?: string;
}

export interface MdesOptions {
    alphaDeg: number;
    symmetric?: boolean;
    filterStrength?: number;
    /** Target velocity distribution (nc circle-plane points). Omit to use current airfoil's. */
    targetQ?: number[];
    /** Number of circle-plane points (must be 2^n + 1, e.g. 129). Default 129. */
    nc?: number;
}

/**
 * Run full-inverse design (MDES) via circle-plane conformal mapping.
 *
 * This is XFOIL's full-inverse method. It maps the airfoil to a circle plane
 * via Fourier coefficients, then recovers geometry from a target Q(w).
 */
export function runFullInverseDesign(
    coordinates: { x: number; y: number }[],
    options: MdesOptions,
): MdesDesignResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const coordsFlat = pointsToFlat(coordinates);
    const nc = options.nc ?? 129;
    const targetQ = options.targetQ ? new Float64Array(options.targetQ) : new Float64Array(0);

    return full_inverse_design_mdes(
        coordsFlat,
        options.alphaDeg,
        options.symmetric ?? false,
        options.filterStrength ?? 0,
        targetQ,
        nc,
    ) as MdesDesignResult;
}

// ============================================================================
// Multi-Element Analysis
// ============================================================================

export interface MultiElementBodyInput {
    name: string;
    coords: { x: number; y: number }[];
    position?: { x: number; y: number; angle: number };
}

export interface PerBodyResult {
    name: string;
    cl: number;
    cp: number[];
    cp_x: number[];
    gamma: number[];
}

export interface MultiElementResult {
    cl_total: number;
    cm_total: number;
    per_body: PerBodyResult[];
    success: boolean;
    error?: string;
}

/**
 * Run multi-element inviscid panel method analysis.
 *
 * Each body is defined by its coordinates and an optional position transform
 * (translation + rotation about leading edge).
 */
export function analyzeMultiElement(
    bodies: MultiElementBodyInput[],
    alphaDeg: number,
): MultiElementResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const bodiesJson = JSON.stringify(bodies.map(b => ({
        name: b.name,
        coords: b.coords.flatMap(p => [p.x, p.y]),
        position: b.position ?? { x: 0, y: 0, angle: 0 },
    })));

    return analyze_multi_element_wasm(bodiesJson, alphaDeg) as MultiElementResult;
}

export interface PerBodyViscousResult {
    name: string;
    cl: number;
    cd: number;
    cm: number;
    converged: boolean;
    success: boolean;
    error?: string;
}

export interface MultiElementViscousResult {
    cl_total: number;
    cd_total: number;
    cm_total: number;
    per_body: PerBodyViscousResult[];
    all_converged: boolean;
    success: boolean;
    error?: string;
}

/**
 * Run multi-element viscous analysis (Phase 1: independent per-body BL).
 *
 * Each body runs its own viscous solver independently. No cross-body
 * viscous coupling. Suitable for well-separated elements.
 */
export function analyzeMultiElementViscous(
    bodies: MultiElementBodyInput[],
    alphaDeg: number,
    reynolds: number = 1e6,
    mach: number = 0,
    ncrit: number = 9,
    maxIterations: number = 100,
): MultiElementViscousResult {
    if (!initialized) {
        throw new Error('WASM not initialized. Call initWasm() first.');
    }

    const bodiesJson = JSON.stringify(bodies.map(b => ({
        name: b.name,
        coords: b.coords.flatMap(p => [p.x, p.y]),
        position: b.position ?? { x: 0, y: 0, angle: 0 },
    })));

    return analyze_multi_element_viscous_wasm(
        bodiesJson, alphaDeg, reynolds, mach, ncrit, maxIterations
    ) as MultiElementViscousResult;
}

// Re-export the RustFoil class and WasmSmokeSystem for advanced usage
export { RustFoil, WasmSmokeSystem };
