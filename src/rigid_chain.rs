//! Generic N-DOF pure rigid serial chain (every joint is a revolute screw joint, jointType==1
//! only, no flexible beam, no PD control). Used for the optimal-control examples (planar 2-link
//! pendulum, rigid 3-link lab robot), where the control input `u` is a free decision variable
//! rather than a PD law. Deliberately kept separate from kinematics.rs/dynamics.rs (the
//! 16-DOF flexible-beam simulation, already validated against MATLAB): those two systems have
//! different force terms (PD + beam stiffness vs. free control), and keeping them apart means
//! this new code can't regress what already works. Both reuse the same se3.rs primitives.

use crate::se3::{self, Vec6};
use nalgebra::{Matrix4, SMatrix, SVector, Vector3};

pub struct ChainDef<const N: usize> {
    pub screw_axis: [Vec6; N],
    pub g_ref: [Matrix4<f64>; N],
    pub parent: [i32; N], // -1 = world
    pub m_gen_diag: [[f64; 6]; N],
    pub mass: [f64; N],
    pub g_const: f64,
}

pub type Jac<const N: usize> = SMatrix<f64, 6, N>;
pub type QVec<const N: usize> = SVector<f64, N>;

fn mgen6<const N: usize>(chain: &ChainDef<N>, i: usize) -> nalgebra::Vector6<f64> {
    nalgebra::Vector6::from_row_slice(&chain.m_gen_diag[i])
}

pub fn g_rel<const N: usize>(chain: &ChainDef<N>, i: usize, q: &[f64; N]) -> Matrix4<f64> {
    chain.g_ref[i] * se3::exp_se3_screw(&chain.screw_axis[i], q[i])
}

pub fn all_g_rel<const N: usize>(chain: &ChainDef<N>, q: &[f64; N]) -> [Matrix4<f64>; N] {
    core::array::from_fn(|i| g_rel(chain, i, q))
}

pub fn fk<const N: usize>(chain: &ChainDef<N>, g_rel_all: &[Matrix4<f64>; N]) -> [Matrix4<f64>; N] {
    let mut g = [Matrix4::identity(); N];
    for i in 0..N {
        let parent_g = if chain.parent[i] < 0 {
            Matrix4::identity()
        } else {
            g[chain.parent[i] as usize]
        };
        g[i] = parent_g * g_rel_all[i];
    }
    g
}

pub fn jacobian<const N: usize>(
    chain: &ChainDef<N>,
    g_rel_all: &[Matrix4<f64>; N],
) -> [Jac<N>; N] {
    let mut j: [Jac<N>; N] = [Jac::<N>::zeros(); N];
    for i in 0..N {
        let mut ji = if chain.parent[i] < 0 {
            Jac::<N>::zeros()
        } else {
            se3::ad_inv(&g_rel_all[i]) * j[chain.parent[i] as usize]
        };
        ji.set_column(i, &chain.screw_axis[i]);
        j[i] = ji;
    }
    j
}

fn discrete_velocities<const N: usize>(
    chain: &ChainDef<N>,
    g_a: &[Matrix4<f64>; N],
    g_b: &[Matrix4<f64>; N],
    h: f64,
) -> [Vec6; N] {
    let mut eta = [Vec6::zeros(); N];
    for i in 0..N {
        let term = if chain.parent[i] < 0 {
            se3::se3_inv(&g_a[i]) * g_b[i]
        } else {
            se3::se3_inv(&g_a[i]) * se3::cay_se3(&(h * eta[chain.parent[i] as usize])) * g_b[i]
        };
        eta[i] = se3::cay_inv_se3(&term) / h;
    }
    eta
}

fn weighted_momentum<const N: usize>(chain: &ChainDef<N>, eta: &Vec6, h: f64, frame: usize, sign: f64) -> Vec6 {
    let t = se3::cay_rtd_inv_se3(&(*eta * (sign * h)));
    let mv = eta.component_mul(&mgen6(chain, frame));
    t.transpose() * mv
}

fn gravity_wrench<const N: usize>(chain: &ChainDef<N>, g_abs: &[Matrix4<f64>; N]) -> [Vec6; N] {
    let up = Vector3::new(0.0, 0.0, 1.0);
    let mut f = [Vec6::zeros(); N];
    for i in 0..N {
        let r = se3::rot(&g_abs[i]);
        let lin = chain.g_const * chain.mass[i] * (r.transpose() * up);
        f[i] = se3::vec6(Vector3::zeros(), lin);
    }
    f
}

struct StepFrozen<const N: usize> {
    g_rel_k: [Matrix4<f64>; N],
    j: [Jac<N>; N],
    f_gen: QVec<N>,
    f_frame: [Vec6; N],
    q_k: [f64; N],
}

fn residual<const N: usize>(chain: &ChainDef<N>, frozen: &StepFrozen<N>, q_k1: &[f64; N], h: f64) -> QVec<N> {
    let g_rel_k1 = all_g_rel(chain, q_k1);
    let eta_k = discrete_velocities(chain, &frozen.g_rel_k, &g_rel_k1, h);
    let mut res = frozen.f_gen;
    for i in 0..N {
        let own = weighted_momentum(chain, &eta_k[i], h, i, 1.0) + frozen.f_frame[i];
        res += frozen.j[i].transpose() * own;
    }
    res
}

fn compute_m_beam<const N: usize>(chain: &ChainDef<N>, j: &[Jac<N>; N], eta_k: &[Vec6; N], h: f64) -> SMatrix<f64, N, N> {
    let mut m = SMatrix::<f64, N, N>::zeros();
    for i in 0..N {
        let a = se3::cay_rtd_inv_se3(&(eta_k[i] * h)).transpose();
        let d = SMatrix::<f64, 6, 6>::from_diagonal(&mgen6(chain, i));
        m += j[i].transpose() * a * d * j[i];
    }
    m
}

fn to_arr<const N: usize>(v: &QVec<N>) -> [f64; N] {
    core::array::from_fn(|i| v[i])
}
fn to_vec<const N: usize>(a: &[f64; N]) -> QVec<N> {
    QVec::<N>::from_row_slice(a)
}

const TOL: f64 = 1e-11;
const MAX_ITER: usize = 100;

fn solve_del<const N: usize>(
    chain: &ChainDef<N>,
    frozen: &StepFrozen<N>,
    mut q_guess: [f64; N],
    h: f64,
    force_iterate: bool,
) -> [f64; N] {
    let mut res = residual(chain, frozen, &q_guess, h);
    if res.norm() <= TOL && !force_iterate {
        return q_guess;
    }

    let eta_k = {
        let g_rel_k1 = all_g_rel(chain, &q_guess);
        discrete_velocities(chain, &frozen.g_rel_k, &g_rel_k1, h)
    };
    let m_beam = compute_m_beam(chain, &frozen.j, &eta_k, h);
    let mut h_inv = (m_beam / h)
        .try_inverse()
        .expect("DEL Jacobian approx singular (degenerate rigid-chain configuration)");

    let mut iters = 0usize;
    loop {
        let q_old = to_vec(&q_guess);
        let res_old = res;
        let step = h_inv * res;
        q_guess = to_arr(&(q_old - step));
        iters += 1;

        res = residual(chain, frozen, &q_guess, h);
        if res.norm() <= TOL || iters >= MAX_ITER {
            break;
        }

        let s = to_vec(&q_guess) - q_old;
        let y = res - res_old;
        let hy = h_inv * y;
        let denom = s.dot(&hy);
        if denom.abs() > 1e-300 {
            let num = (s - hy) * (s.transpose() * h_inv);
            h_inv += num / denom;
        }
    }
    q_guess
}

/// Forward-simulate the chain under free control inputs `u` (one `[f64;N]` per step, B = identity
/// selection on every DOF, no PD/damping), starting at rest from `q0` (qDot0 assumed zero, as in
/// every ELARA optimal-control example). Returns the full `q` trajectory, length `u.len()+1`.
pub fn simulate<const N: usize>(chain: &ChainDef<N>, h: f64, q0: [f64; N], u: &[[f64; N]]) -> Vec<[f64; N]> {
    let mut traj = Vec::with_capacity(u.len() + 1);
    traj.push(q0);
    if u.is_empty() {
        return traj;
    }

    let g_rel_0 = all_g_rel(chain, &q0);
    let g_0 = fk(chain, &g_rel_0);
    let j_0 = jacobian(chain, &g_rel_0);
    let grav_0 = gravity_wrench(chain, &g_0);

    let a0 = 0.5;
    let mut f_gen_0 = QVec::<N>::zeros();
    for i in 0..N {
        f_gen_0[i] = -h * (1.0 - a0) * u[0][i];
    }
    let mut f_frame_0 = [Vec6::zeros(); N];
    for i in 0..N {
        f_frame_0[i] = grav_0[i] * (h * (1.0 - a0));
    }
    let frozen_0 = StepFrozen { g_rel_k: g_rel_0, j: j_0, f_gen: f_gen_0, f_frame: f_frame_0, q_k: q0 };
    let q1 = solve_del(chain, &frozen_0, q0, h, true);
    traj.push(q1);

    let mut q_k0 = q0;
    let mut q_k = q1;
    let mut g_rel_k0 = g_rel_0;

    for u_k in u.iter().skip(1) {
        let g_rel_k = all_g_rel(chain, &q_k);
        let g_k = fk(chain, &g_rel_k);
        let j_k = jacobian(chain, &g_rel_k);
        let eta_k0 = discrete_velocities(chain, &g_rel_k0, &g_rel_k, h);
        let grav_k = gravity_wrench(chain, &g_k);

        let mut f_gen = QVec::<N>::zeros();
        for i in 0..N {
            f_gen[i] = -h * u_k[i];
        }
        let mut f_frame = [Vec6::zeros(); N];
        for i in 0..N {
            f_frame[i] = grav_k[i] * h - weighted_momentum(chain, &eta_k0[i], h, i, -1.0);
        }

        let frozen = StepFrozen { g_rel_k, j: j_k, f_gen, f_frame, q_k };
        let mut guess = [0.0; N];
        for i in 0..N {
            guess[i] = 2.0 * q_k[i] - q_k0[i];
        }
        let q_k1 = solve_del(chain, &frozen, guess, h, false);
        traj.push(q_k1);

        g_rel_k0 = g_rel_k;
        q_k0 = q_k;
        q_k = q_k1;
    }

    traj
}
