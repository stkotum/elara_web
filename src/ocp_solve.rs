//! Direct-shooting trajectory optimizer for a `rigid_chain::ChainDef<N>`: free controls
//! `u_0..u_{M-1}` are the only decision variables (states are eliminated via forward simulation,
//! not treated as separate collocation variables the way ELARA's CasADi/IPOPT transcription
//! does - see reference/ocp_math_spec.md for why that's a faithful reformulation of the same
//! underlying optimal-control problem, not a different one). The running cost here
//! (`w_u * h * sum(||u_k||^2)`) is exactly quadratic in `u`, so each SQP iteration solves an
//! exact QP subproblem (finite-difference-linearized nonlinear boundary constraints, exact
//! quadratic cost) via the Schur complement on the small (2N x 2N) constraint-multiplier system,
//! not a generic sparse KKT solve.

use crate::rigid_chain::{self, ChainDef};
use nalgebra::{DMatrix, DVector};

pub struct OcpProblem<const N: usize> {
    pub h: f64,
    pub q0: [f64; N],
    pub q_final_target: [f64; N],
    pub qdot_final_target: [f64; N],
    pub n_steps: usize,
    pub w_u: f64,
}

pub struct OcpResult<const N: usize> {
    pub u: Vec<[f64; N]>,
    pub q_traj: Vec<[f64; N]>,
    pub cost: f64,
    pub constraint_violation: f64,
    pub iterations: usize,
}

fn unflatten<const N: usize>(v: &DVector<f64>, n_steps: usize) -> Vec<[f64; N]> {
    (0..n_steps).map(|k| core::array::from_fn(|i| v[k * N + i])).collect()
}

fn boundary_residual<const N: usize>(chain: &ChainDef<N>, prob: &OcpProblem<N>, u: &[[f64; N]]) -> DVector<f64> {
    let traj = rigid_chain::simulate(chain, prob.h, prob.q0, u);
    let q_final = traj[traj.len() - 1];
    let q_prev = traj[traj.len() - 2];
    let mut g = DVector::zeros(2 * N);
    for i in 0..N {
        g[i] = q_final[i] - prob.q_final_target[i];
        let qdot_final = (q_final[i] - q_prev[i]) / prob.h;
        g[N + i] = qdot_final - prob.qdot_final_target[i];
    }
    g
}

pub fn solve<const N: usize>(chain: &ChainDef<N>, prob: &OcpProblem<N>) -> OcpResult<N> {
    let nvar = prob.n_steps * N;
    let mut u_vec = DVector::<f64>::zeros(nvar);
    let eps = 1e-6;
    let h_diag = 2.0 * prob.w_u * prob.h; // constant Hessian diagonal of the exactly-quadratic cost
    let max_outer = 40;
    let mut iters = 0usize;

    for outer in 0..max_outer {
        iters = outer;
        let u = unflatten::<N>(&u_vec, prob.n_steps);
        let g = boundary_residual(chain, prob, &u);
        let gnorm = g.norm();
        if gnorm < 1e-9 && outer > 0 {
            break;
        }

        let grad_j = &u_vec * h_diag;

        let mcon = 2 * N;
        let mut jg = DMatrix::<f64>::zeros(mcon, nvar);
        for j in 0..nvar {
            let mut up = u_vec.clone();
            up[j] += eps;
            let mut um = u_vec.clone();
            um[j] -= eps;
            let gp = boundary_residual(chain, prob, &unflatten::<N>(&up, prob.n_steps));
            let gm = boundary_residual(chain, prob, &unflatten::<N>(&um, prob.n_steps));
            jg.set_column(j, &((gp - gm) / (2.0 * eps)));
        }

        // Schur complement: H = h_diag * I, so H^-1 is just a scalar divide.
        // KKT: du = -H^-1(grad_J + Jg'*lambda), g + Jg*du = 0
        //   => (Jg H^-1 Jg') lambda = g - Jg H^-1 grad_J
        let jg_hinv = jg.clone() / h_diag;
        let schur = &jg_hinv * jg.transpose();
        let rhs = &g - &jg_hinv * &grad_j;
        let lambda = schur
            .clone()
            .lu()
            .solve(&rhs)
            .expect("OCP KKT Schur system singular (degenerate constraint Jacobian)");

        let du = (-&grad_j - jg.transpose() * &lambda) / h_diag;

        let mut alpha = 1.0;
        for _ in 0..20 {
            let u_try = &u_vec + &du * alpha;
            let g_try = boundary_residual(chain, prob, &unflatten::<N>(&u_try, prob.n_steps));
            if g_try.norm() <= gnorm + 1e-12 || gnorm < 1e-6 {
                u_vec = u_try;
                break;
            }
            alpha *= 0.5;
        }
    }

    let u_final = unflatten::<N>(&u_vec, prob.n_steps);
    let traj = rigid_chain::simulate(chain, prob.h, prob.q0, &u_final);
    let g_final = boundary_residual(chain, prob, &u_final);
    let cost = prob.w_u * prob.h * u_vec.iter().map(|x| x * x).sum::<f64>();

    OcpResult { u: u_final, q_traj: traj, cost, constraint_violation: g_final.norm(), iterations: iters }
}
