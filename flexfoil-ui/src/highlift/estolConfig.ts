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

/**
 * Cove-cutout parameters for the main element (chord-normalized). The cove
 * follows the airfoil's own upper and lower surfaces (offset inward by the
 * blunt-TE thickness) from each cut forward to `coveX`, joined by a vertical back
 * wall at `coveX` with sharp corners. `coveX` sits forward of both cuts.
 */
export interface CoveCutouts {
  /** x where the upper surface is cut (the lip). */
  upperCutX: number;
  /** x where the lower surface is cut. */
  lowerCutX: number;
  /** x of the cove back wall (forward of both cuts). */
  coveX: number;
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
 * Flaperon shape, derived from the main airfoil rather than an independent foil.
 * Aft of `upperCutX` / `lowerCutX` the surfaces are exactly LS(1)-0417 (blunt-
 * truncated at the TE like the other elements). Forward of the cuts the nose is a
 * 5-control-point clamped cubic B-spline whose end tangents match the airfoil
 * slope at each cut: P0 = upper cut, P1 = P0 + noseTangentUpper·tHat(upperCutX),
 * P2 = noseTip, P3 = P4 + noseTangentLower·tHat(lowerCutX), P4 = lower cut.
 * Defined in the airfoil (global) frame; deployment is applied on top.
 */
/**
 * Vane shape: a free-form outline (no longer a NACA foil). The outline is a closed
 * periodic cubic B-spline through 5 control points — the trailing-edge location
 * `te` plus 4 other `points` (loop order [te, ...points]) — with a blunt TE of the
 * airfoil's `bluntThickness` cut at `te` (shared with the cove lips and flaperon).
 */
export interface VaneShape {
  /** Trailing-edge location; the blunt TE is cut here. */
  te: V2;
  /** The 4 other B-spline control points (with `te`, a 5-point closed loop). */
  points: V2[];
}

export interface FlaperonShape {
  /** Upper-surface cut x; aft of it the upper surface = the airfoil. */
  upperCutX: number;
  /** Lower-surface cut x; aft of it the lower surface = the airfoil. */
  lowerCutX: number;
  /** P1 distance from the upper cut along the forward airfoil tangent. */
  noseTangentUpper: number;
  /** P3 distance from the lower cut along the forward airfoil tangent. */
  noseTangentLower: number;
  /** P2, the free middle nose control point (global coordinates). */
  noseTip: V2;
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
  /** Maximum deployment (the deploy slider's upper bound), in chord units. */
  length: number;
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
  /** Vane shape (free-form B-spline). */
  vane: VaneShape;
  /** Flaperon shape, derived from the main airfoil. */
  flaperon: FlaperonShape;
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
      upperCutX: 0.8,
      lowerCutX: 0.68,
      coveX: 0.55,
    },
    vane: {
      te: [0.75, 0.05],
      points: [[0.655, 0.055], [0.555, 0.015], [0.55, -0.08], [0.66, 0.045]],
    },
    flaperon: {
      upperCutX: 0.9,
      lowerCutX: 0.75,
      noseTangentUpper: 0.025,
      noseTangentLower: 0.11,
      noseTip: [0.715, 0.075],
    },
    track: { anchor: [0.675, -0.03], angleDeg: -5.5, linearLength: 0.1, arcRadius: -0.15, length: 0.25 },
    axelSpacing: 0.135,
    flaperonHinge: { pivot: [0.76, -0.02] },
  },
  operation: {
    deploy: 0.0,
    flaperonAngleDeg: 0.0,
  },
};
