//! Cross-checks the from-scratch Rust OCP solver against a REAL CasADi+IPOPT solve of
//! `optimal_control_planar_robot.m` (reference/ocp_planar/, captured via `matlab -batch` with a
//! locally installed CasADi). The Rust port uses direct shooting (see ocp_solve.rs) rather than
//! MATLAB's direct collocation, so exact trajectory agreement isn't expected - but both solve
//! the same continuous optimal-control problem, so the optimal cost should land close.

use elara_core::ocp_solve::{solve, OcpProblem};
use elara_core::systems::planar_robot;

#[test]
fn matches_matlab_ipopt_cost_and_boundary() {
    let chain = planar_robot::def();
    let prob = OcpProblem::<{ planar_robot::N }> {
        h: 0.01,
        q0: [std::f64::consts::FRAC_PI_2, 0.0],
        q_final_target: [-std::f64::consts::FRAC_PI_2, 0.0],
        qdot_final_target: [0.0, 0.0],
        n_steps: 100,
        w_u: 0.5,
    };

    let result = solve(&chain, &prob);

    eprintln!("iterations: {}", result.iterations);
    eprintln!("constraint violation: {:e}", result.constraint_violation);
    eprintln!("cost (Rust shooting): {}", result.cost);
    eprintln!("MATLAB IPOPT cost: 80.5257");
    let final_q = result.q_traj.last().unwrap();
    eprintln!("final q: {:?} (target [-1.5708, 0])", final_q);

    assert!(
        result.constraint_violation < 1e-6,
        "boundary conditions not satisfied: {:e}",
        result.constraint_violation
    );
    assert!(
        (final_q[0] - prob.q_final_target[0]).abs() < 1e-4,
        "q1 final mismatch: {}",
        final_q[0]
    );
    assert!(
        (final_q[1] - prob.q_final_target[1]).abs() < 1e-4,
        "q2 final mismatch: {}",
        final_q[1]
    );

    // Same continuous problem, different discretization/transcription (shooting vs.
    // collocation) - expect the optimal cost to be close, not bit-identical.
    let matlab_cost = 80.5257_f64;
    let rel_err = (result.cost - matlab_cost).abs() / matlab_cost;
    eprintln!("relative cost error vs MATLAB: {:.4}", rel_err);
    assert!(rel_err < 0.1, "cost differs from MATLAB by more than 10%: rel_err={:.4}", rel_err);
}
