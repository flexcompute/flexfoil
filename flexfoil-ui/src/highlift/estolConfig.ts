/**
 * Configuration for the 3-element high-lift airfoil, organized as a 3-level
 * hierarchy that mirrors the design dependency chain:
 *
 *   L1  main      — the base airfoil (fixed LS(1)-0417). Changing it invalidates
 *                   everything downstream.
 *   L2  design    — the high-lift element design built on that airfoil: main cove
 *                   cutout, vane & flaperon shapes, flap track, flaperon hinge.
 *   L3  operation — the operating point on a fixed design: deployment + flaperon angle.
 *
 * Single source of truth for the default geometry consumed by the /highlift
 * design tool and by the geometry engine in ./geometry.ts. Originally ported from
 * ~/himalaya/geometry/airfoils/estol_config.yaml.
 */

export type V2 = [number, number];

/** Explicit upper/lower surface coordinates for a tabulated airfoil (e.g. LS(1)-0417). */
export interface AirfoilCoords {
  /** LE → TE, ascending x. */
  upper: V2[];
  /** LE → TE, ascending x. */
  lower: V2[];
}

/** Cove-cutout parameters for the main element (chord-normalized). */
export interface CoveCutouts {
  /** x where the lower surface is cut to start the cove. */
  lowerCutX: number;
  /** x where the upper surface lip is cut. */
  upperLipCutX: number;
  /** Cove vertex (deepest interior corner) x. */
  coveVertexX: number;
  /** Cove vertex y. */
  coveVertexY: number;
  /** Fillet radius applied at the cove vertex. */
  coveFilletRadius: number;
}

/** A NACA 4-digit element anchored between an LE and TE point. */
export interface NacaElementCfg {
  /** 4-digit NACA designation, e.g. "9621". */
  naca: string;
  /** Stowed leading-edge anchor (global coordinates). */
  stowedLe: V2;
  /** Stowed trailing-edge anchor (global coordinates). */
  stowedTe: V2;
}

/**
 * Flap track: a linear segment from `anchor` along `angleDeg`, followed by a
 * circular arc tangent to it at the junction. The vane+flaperon assembly rides the
 * track on two axels spaced `axelSpacing` apart (see geometry.ts), giving a
 * 1-DOF deployment driven by a single slider.
 */
export interface TrackConfig {
  /** Start of the linear segment (global coords). */
  anchor: V2;
  /** Linear segment direction, degrees. */
  angleDeg: number;
  /** Length of the linear segment before the arc begins. */
  linearLength: number;
  /** Signed arc radius (negative curves the track / assembly TE-down). */
  arcRadius: number;
}

// ---------------------------------------------------------------------------
// Level 1 — main airfoil (foundational; changing it invalidates everything below)
// ---------------------------------------------------------------------------

export interface MainAirfoil {
  /** Base airfoil coordinates (fixed to LS(1)-0417 for now). */
  coords: AirfoilCoords;
  /** Flat blunt TE thickness, in overall-chord units, applied to all elements. */
  bluntThickness: number;
}

// ---------------------------------------------------------------------------
// Level 2 — high-lift element design (the system geometry on the main airfoil)
// ---------------------------------------------------------------------------

/**
 * Flaperon hinge geometry: the pivot (at stowed deployment, global coords) about
 * which the flaperon rotates relative to the vane+flaperon assembly. The rotation
 * angle itself is an operating-point quantity and lives in {@link Operation}.
 */
export interface FlaperonHinge {
  /** Pivot at stowed deployment (global coords); carried by the track motion. */
  pivot: V2;
}

export interface HighLiftDesign {
  /** Main element cove cutout. */
  cove: CoveCutouts;
  /** Vane element shape + stowed anchors. */
  vane: NacaElementCfg;
  /** Flaperon element shape + stowed anchors. */
  flaperon: NacaElementCfg;
  /** Flap track shared by the vane + flaperon assembly. */
  track: TrackConfig;
  /** Chord distance between the two track axels. */
  axelSpacing: number;
  /** Flaperon hinge pivot (geometry only; the angle is in operation). */
  flaperonHinge: FlaperonHinge;
}

// ---------------------------------------------------------------------------
// Level 3 — operating point on a fixed design
// ---------------------------------------------------------------------------

export interface Operation {
  /** Deployment: travel beyond stowed along the track; sLead = design.axelSpacing + deploy. */
  deploy: number;
  /** Flaperon rotation about its hinge, relative to the assembly (deg, negative = TE down). */
  flaperonAngleDeg: number;
}

export interface HighLiftAirfoil {
  main: MainAirfoil;
  design: HighLiftDesign;
  operation: Operation;
}

/** NASA LS(1)-0417 ("GA(W)-1") tabulated coordinates. */
export const LS1_0417: AirfoilCoords = {
  upper: [
    [0.0, 0.0], [0.005, 0.0245], [0.0125, 0.0366], [0.025, 0.05], [0.05, 0.0678],
    [0.1, 0.0909], [0.15, 0.1065], [0.2, 0.1165], [0.25, 0.1225], [0.3, 0.1255],
    [0.35, 0.126], [0.4, 0.1245], [0.45, 0.1215], [0.5, 0.1165], [0.55, 0.1095],
    [0.6, 0.101], [0.65, 0.091], [0.7, 0.0795], [0.75, 0.0665], [0.8, 0.052],
    [0.85, 0.0365], [0.9, 0.0205], [0.95, 0.0085], [0.98, 0.0025], [1.0, 0.0],
  ],
  lower: [
    [0.0, 0.0], [0.005, -0.0152], [0.0125, -0.025], [0.025, -0.0345], [0.05, -0.0465],
    [0.1, -0.0585], [0.15, -0.0645], [0.2, -0.0675], [0.25, -0.0685], [0.3, -0.0685],
    [0.35, -0.0675], [0.4, -0.0655], [0.45, -0.0625], [0.5, -0.058], [0.55, -0.052],
    [0.6, -0.044], [0.65, -0.034], [0.7, -0.021], [0.75, -0.011], [0.8, -0.004],
    [0.85, 0.0], [0.9, 0.002], [0.95, 0.0015], [0.98, 0.0005], [1.0, 0.0],
  ],
};

/** Default configuration matching the reference paper geometry. */
export const DEFAULT_HIGH_LIFT_AIRFOIL: HighLiftAirfoil = {
  main: {
    coords: LS1_0417,
    bluntThickness: 0.003,
  },
  design: {
    cove: {
      lowerCutX: 0.51,
      upperLipCutX: 0.85,
      coveVertexX: 0.62,
      coveVertexY: 0.055,
      coveFilletRadius: 0.2,
    },
    vane: { naca: '9621', stowedLe: [0.55, -0.048], stowedTe: [0.73, 0.035] },
    flaperon: { naca: '6311', stowedLe: [0.69, -0.011], stowedTe: [1.0, 0.0] },
    track: { anchor: [0.675, -0.03], angleDeg: -6.0, linearLength: 0.155, arcRadius: -0.135 },
    axelSpacing: 0.095,
    flaperonHinge: { pivot: [0.69, -0.011] },
  },
  operation: {
    deploy: 0.12,
    flaperonAngleDeg: 0.0,
  },
};
