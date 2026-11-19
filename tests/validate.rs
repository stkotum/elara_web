//! Validates the Rust port against a real MATLAB run of
//! `simulation_rigid_flexible_robot_PD_control.m` (reference/vi_q_qdot.csv), generated via
//! `matlab -batch` from the actual installed ELARA toolbox. See reference/math_spec.md.

use elara_core::dynamics::{simulate, SolverConfig};
use elara_core::system_data::N_DOF;

fn load_reference(path: &str) -> Vec<(f64, [f64; N_DOF], [f64; N_DOF])> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .expect("reference CSV missing");
    let mut rows = Vec::new();
    for rec in rdr.records() {
        let rec = rec.unwrap();
        let vals: Vec<f64> = rec.iter().map(|s| s.parse().unwrap()).collect();
        assert_eq!(vals.len(), 1 + 2 * N_DOF);
        let t = vals[0];
        let mut q = [0.0; N_DOF];
        let mut qd = [0.0; N_DOF];
        q.copy_from_slice(&vals[1..1 + N_DOF]);
        qd.copy_from_slice(&vals[1 + N_DOF..1 + 2 * N_DOF]);
        rows.push((t, q, qd));
    }
    rows
}

#[test]
fn matches_matlab_variational_integrator() {
    let reference = load_reference("reference/vi_q_qdot.csv");
    let t_end = reference.last().unwrap().0;

    let cfg = SolverConfig::default();
    let traj = simulate(t_end, [0.0; N_DOF], [0.0; N_DOF], &cfg);

    assert_eq!(
        traj.q.len(),
        reference.len(),
        "step count mismatch: port={} matlab={}",
        traj.q.len(),
        reference.len()
    );

    let mut max_abs_err = 0.0f64;
    let mut max_abs_err_at = (0usize, 0usize);
    let mut final_err = [0.0f64; N_DOF];

    for (k, (t_ref, q_ref, _qd_ref)) in reference.iter().enumerate() {
        assert!((traj.t[k] - t_ref).abs() < 1e-9, "time grid mismatch at step {k}");
        for d in 0..N_DOF {
            let e = (traj.q[k][d] - q_ref[d]).abs();
            if e > max_abs_err {
                max_abs_err = e;
                max_abs_err_at = (k, d);
            }
            if k == reference.len() - 1 {
                final_err[d] = e;
            }
        }
    }

    eprintln!("max |dq| = {max_abs_err:e} at step {} dof {}", max_abs_err_at.0, max_abs_err_at.1);
    eprintln!("final-step error per DOF: {final_err:?}");
    eprintln!("final q (port):   {:?}", traj.q.last().unwrap());
    eprintln!("final q (matlab): {:?}", reference.last().unwrap().1);

    // Observed max |dq| over the full 5s/1280-step/16-DOF trajectory is ~5e-10 (float-order
    // roundoff between MATLAB and Rust, well inside the solver's own 5e-12 residual tolerance
    // compounding over steps). 1e-6 leaves ample margin while still catching real regressions.
    assert!(
        max_abs_err < 1e-6,
        "port diverges from MATLAB reference: max |dq| = {max_abs_err:e}"
    );
}
