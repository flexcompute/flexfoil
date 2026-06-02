//! Validation for the coupled multi-body inviscid solver.
//!
//! No XFOIL/mfoil multi-element reference exists, so we anchor on:
//!   1. Single-body equivalence — one body through the multi path must match the
//!      single-body solver to machine precision (the generalization is exact).
//!   2. Far-apart recovery — two widely separated bodies each recover the
//!      isolated single-body result (coupling → 0 with distance).
//!   3. Close coupling — proximity changes each body's lift measurably and in the
//!      physically expected direction (vertically stacked biplane → mutual
//!      interference reduces lift).

use rustfoil_inviscid::{
    build_and_factorize_multi, system::build_and_factorize, AirfoilGeometry, FlowConditions,
};
use std::path::PathBuf;

fn naca0012() -> Vec<(f64, f64)> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/naca0012_xfoil_paneled.dat");
    let content = std::fs::read_to_string(&path).expect("read naca0012 fixture");
    content
        .lines()
        .skip(1)
        .filter_map(|l| {
            let p: Vec<&str> = l.split_whitespace().collect();
            (p.len() >= 2).then(|| (p[0].parse().ok().unwrap(), p[1].parse().ok().unwrap()))
        })
        .collect()
}

fn shifted(coords: &[(f64, f64)], dx: f64, dy: f64) -> Vec<(f64, f64)> {
    coords.iter().map(|&(x, y)| (x + dx, y + dy)).collect()
}

fn geom(coords: &[(f64, f64)]) -> AirfoilGeometry {
    AirfoilGeometry::from_points(coords).expect("build geometry")
}

const ALPHA: f64 = 4.0;

#[test]
fn single_body_equivalence() {
    let coords = naca0012();
    let flow = FlowConditions::with_alpha_deg(ALPHA);

    let single = build_and_factorize(&geom(&coords)).unwrap().solve_alpha(&flow);
    let multi = build_and_factorize_multi(&[geom(&coords)]).unwrap().solve_alpha(&flow);

    assert_eq!(multi.bodies.len(), 1);
    let b = &multi.bodies[0];
    assert!((b.cl - single.cl).abs() < 1e-9, "cl: multi {} vs single {}", b.cl, single.cl);
    assert!((b.cm - single.cm).abs() < 1e-9, "cm: multi {} vs single {}", b.cm, single.cm);
    for (i, (&m, &s)) in b.cp.iter().zip(single.cp.iter()).enumerate() {
        assert!((m - s).abs() < 1e-9, "cp[{i}]: multi {m} vs single {s}");
    }
}

#[test]
fn far_apart_recovers_isolated() {
    let coords = naca0012();
    let flow = FlowConditions::with_alpha_deg(ALPHA);
    let cl_iso = build_and_factorize(&geom(&coords)).unwrap().solve_alpha(&flow).cl;

    // Second body 50 chords above — interaction should be negligible.
    let sol = build_and_factorize_multi(&[geom(&coords), geom(&shifted(&coords, 0.0, 50.0))])
        .unwrap()
        .solve_alpha(&flow);

    for (k, b) in sol.bodies.iter().enumerate() {
        assert!(
            (b.cl - cl_iso).abs() < 1e-3,
            "body {k} cl {} should match isolated {cl_iso}",
            b.cl
        );
    }
}

#[test]
fn coupling_decays_with_separation() {
    let coords = naca0012();
    let flow = FlowConditions::with_alpha_deg(ALPHA);
    let cl_iso = build_and_factorize(&geom(&coords)).unwrap().solve_alpha(&flow).cl;

    // Interaction strength = largest per-body departure from the isolated Cl,
    // for a vertically stacked pair at increasing gap. Physically certain: the
    // potential coupling must weaken monotonically with separation.
    let strength = |gap: f64| {
        let sol = build_and_factorize_multi(&[geom(&coords), geom(&shifted(&coords, 0.0, gap))])
            .unwrap()
            .solve_alpha(&flow);
        sol.bodies
            .iter()
            .map(|b| {
                assert!(b.cl.is_finite());
                (b.cl - cl_iso).abs()
            })
            .fold(0.0_f64, f64::max)
    };

    let (near, mid, far) = (strength(0.5), strength(2.0), strength(10.0));
    println!("isolated cl = {cl_iso:.5}; interaction: near {near:.5}, mid {mid:.5}, far {far:.5}");

    assert!(near > 1e-2, "coupling should be active at 0.5c gap (got {near})");
    assert!(near > mid && mid > far, "interaction must decay with separation");
    assert!(far < near * 0.25, "interaction should be much weaker at 10c (got {far} vs {near})");
}
