use crate::kinematics::{self, Jac};
use crate::se3::{self, Vec6};
use crate::system_data::*;
use nalgebra::{Matrix4, SMatrix, SVector};

type QVec = SVector<f64, N_DOF>;
type MMat = SMatrix<f64, N_DOF, N_DOF>;

fn to_arr(v: &QVec) -> [f64; N_DOF] {
    core::array::from_fn(|i| v[i])
}
fn to_vec(a: &[f64; N_DOF]) -> QVec {
    QVec::from_row_slice(a)
}

fn mgen_diag6(i: usize) -> nalgebra::Vector6<f64> {
    nalgebra::Vector6::from_row_slice(&M_GEN_DIAG[i])
}

/// The 16-DOF generalized stiffness/damping/reference arrays. Entries 0..4 are the "PD gains"
/// on the 4 rigid-ish screw joints (user-adjustable in the interactive demo); entries 4..16 are
/// the flexible beam's own physical stiffness/damping, fixed by its material and geometry.
#[derive(Clone, Copy)]
pub struct Gains {
    pub c_sys: [f64; N_DOF],
    pub d_sys: [f64; N_DOF],
    pub q_ref: [f64; N_DOF],
}

impl Gains {
    /// `cSys[0:4]`, `dSys[0:4]`, `qRef[0:4]` from the PD panel; beam DOF (4..16) always use the
    /// physical values baked into `system_data` (ELARA's own `cSys = l*Ba'*diag(Cgen)` etc.).
    pub fn new(p_gain: [f64; 4], d_gain: [f64; 4], q_des: [f64; 4]) -> Self {
        let mut c_sys = C_SYS;
        let mut d_sys = D_SYS;
        let mut q_ref = Q_REF;
        c_sys[0..4].copy_from_slice(&p_gain);
        d_sys[0..4].copy_from_slice(&d_gain);
        q_ref[0..4].copy_from_slice(&q_des);
        Self { c_sys, d_sys, q_ref }
    }

    /// Exactly the example script's values: P=100, D=30, qDes = [45,-45,90,-80] deg.
    pub fn example_default() -> Self {
        Self::new(
            [100.0; 4],
            [30.0; 4],
            [
                45f64.to_radians(),
                (-45f64).to_radians(),
                90f64.to_radians(),
                (-80f64).to_radians(),
            ],
        )
    }
}

/// Recursive discrete body-fixed velocity, `eta(i)` for the interval between relative-transform
/// sets `g_a` (start) and `g_b` (end), per `computeDiscreteAbsoluteVelocities`.
fn discrete_velocities(
    g_a: &[Matrix4<f64>; N_FRAMES],
    g_b: &[Matrix4<f64>; N_FRAMES],
    h: f64,
) -> [Vec6; N_FRAMES] {
    let mut eta = [Vec6::zeros(); N_FRAMES];
    for i in 0..N_FRAMES {
        let term = if PARENT[i] < 0 {
            se3::se3_inv(&g_a[i]) * g_b[i]
        } else {
            se3::se3_inv(&g_a[i]) * se3::cay_se3(&(h * eta[PARENT[i] as usize])) * g_b[i]
        };
        eta[i] = se3::cay_inv_se3(&term) / h;
    }
    eta
}

/// `cayRTDInvSE3(sign*eta*h)' * MGen * eta`, the recurring "weighted discrete momentum" term.
fn weighted_momentum(eta: &Vec6, h: f64, frame: usize, sign: f64) -> Vec6 {
    let t = se3::cay_rtd_inv_se3(&(*eta * (sign * h)));
    let mv = eta.component_mul(&mgen_diag6(frame));
    t.transpose() * mv
}

/// Gravity-only body-fixed wrench per frame (torque term vanishes: x_a == 0 for every frame in
/// this example, so the point-mass moment-arm cross product is identically zero).
fn gravity_wrench(g_abs: &[Matrix4<f64>; N_FRAMES]) -> [Vec6; N_FRAMES] {
    let up = kinematics::gravity_dir();
    let mut f = [Vec6::zeros(); N_FRAMES];
    for i in 0..N_FRAMES {
        let r = se3::rot(&g_abs[i]);
        let lin = G_CONST * MASS[i] * (r.transpose() * up);
        f[i] = se3::vec6(nalgebra::Vector3::zeros(), lin);
    }
    f
}

/// Everything about a step that is frozen before the Broyden iterations begin (evaluated once
/// at `q_k`, reused across every residual evaluation for this step).
struct StepFrozen {
    g_rel_k: [Matrix4<f64>; N_FRAMES],
    j: [Jac; N_FRAMES],
    f_gen: QVec,
    f_frame: [Vec6; N_FRAMES],
    damping_coeff: f64, // (1 - a)
    q_k: [f64; N_DOF],
}

fn residual(frozen: &StepFrozen, q_k1: &[f64; N_DOF], gains: &Gains, h: f64) -> (QVec, [Vec6; N_FRAMES]) {
    let g_rel_k1 = kinematics::all_g_rel(q_k1);
    let eta_k = discrete_velocities(&frozen.g_rel_k, &g_rel_k1, h);
    let mut res = frozen.f_gen;
    for i in 0..N_DOF {
        res[i] += frozen.damping_coeff * gains.d_sys[i] * (q_k1[i] - frozen.q_k[i]);
    }
    for i in 0..N_FRAMES {
        let own = weighted_momentum(&eta_k[i], h, i, 1.0) + frozen.f_frame[i];
        res += frozen.j[i].transpose() * own;
    }
    (res, eta_k)
}

fn compute_m_beam(j: &[Jac; N_FRAMES], eta_k: &[Vec6; N_FRAMES], h: f64) -> MMat {
    let mut m = MMat::zeros();
    for i in 0..N_FRAMES {
        let a = se3::cay_rtd_inv_se3(&(eta_k[i] * h)).transpose();
        let d = SMatrix::<f64, 6, 6>::from_diagonal(&mgen_diag6(i));
        m += j[i].transpose() * a * d * j[i];
    }
    m
}

fn d_sys_diag(gains: &Gains) -> MMat {
    MMat::from_diagonal(&QVec::from_row_slice(&gains.d_sys))
}

#[derive(Clone, Copy)]
pub struct SolverConfig {
    pub h: f64,
    pub tol: f64,
    pub jac_iter_threshold: usize,
    pub max_iter: usize,
}

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            h: 2f64.powi(-8),
            tol: 5e-12,
            jac_iter_threshold: 5,
            max_iter: 100,
        }
    }
}

/// Broyden "good" quasi-Newton solve of the DEL residual for `q_k1`, carrying the inverse-
/// Jacobian approximation `h_inv` across steps as ELARA's `solveImplicitDELEquBroyden.m` does.
fn solve_del(
    frozen: &StepFrozen,
    mut q_guess: [f64; N_DOF],
    gains: &Gains,
    h_inv: &mut MMat,
    must_recompute: bool,
    force_iterate: bool,
    cfg: &SolverConfig,
) -> ([f64; N_DOF], usize) {
    let (mut res, eta_k) = residual(frozen, &q_guess, gains, cfg.h);

    if res.norm() <= cfg.tol && !force_iterate {
        return (q_guess, 0);
    }

    if must_recompute {
        let m_beam = compute_m_beam(&frozen.j, &eta_k, cfg.h);
        *h_inv = (m_beam / cfg.h + d_sys_diag(gains))
            .try_inverse()
            .expect("DEL Jacobian approx singular");
    }

    let mut iters = 0usize;
    loop {
        let q_old = to_vec(&q_guess);
        let res_old = res;
        let step = *h_inv * res;
        let q_new = q_old - step;
        q_guess = to_arr(&q_new);
        iters += 1;

        let (res_new, _eta_new) = residual(frozen, &q_guess, gains, cfg.h);
        res = res_new;

        if res.norm() <= cfg.tol || iters >= cfg.max_iter {
            break;
        }

        // Broyden "good" rank-1 inverse-Jacobian update.
        let s = to_vec(&q_guess) - q_old;
        let y = res - res_old;
        let hy = *h_inv * y;
        let denom = s.dot(&hy);
        if denom.abs() > 1e-300 {
            let num = (s - hy) * (s.transpose() * *h_inv);
            *h_inv += num / denom;
        }
    }
    (q_guess, iters)
}

fn build_frozen_general(q_k0: &[f64; N_DOF], q_k: &[f64; N_DOF], g_rel_k0: &[Matrix4<f64>; N_FRAMES], gains: &Gains, h: f64) -> (StepFrozen, [Matrix4<f64>; N_FRAMES]) {
    let a_trapez = 0.0; // example setting: fully one-sided damping, see reference/math_spec.md Sec.4
    let g_rel_k = kinematics::all_g_rel(q_k);
    let g_k = kinematics::fk(&g_rel_k);
    let j_k = kinematics::jacobian(q_k, &g_rel_k);
    let eta_k0 = discrete_velocities(g_rel_k0, &g_rel_k, h);
    let grav_k = gravity_wrench(&g_k);

    let mut f_gen = QVec::zeros();
    for i in 0..N_DOF {
        f_gen[i] =
            h * gains.c_sys[i] * (q_k[i] - gains.q_ref[i]) + a_trapez * gains.d_sys[i] * (q_k[i] - q_k0[i]);
    }
    let mut f_frame = [Vec6::zeros(); N_FRAMES];
    for i in 0..N_FRAMES {
        f_frame[i] = grav_k[i] * h - weighted_momentum(&eta_k0[i], h, i, -1.0);
    }

    (
        StepFrozen {
            g_rel_k,
            j: j_k,
            f_gen,
            f_frame,
            damping_coeff: 1.0 - a_trapez,
            q_k: *q_k,
        },
        g_rel_k,
    )
}

/// Live, steppable simulation: holds the two-step state the variational integrator needs
/// between calls, so a UI can call `step()` once per physics tick and re-render.
pub struct SimState {
    pub cfg: SolverConfig,
    pub gains: Gains,
    q_k0: [f64; N_DOF],
    q_k: [f64; N_DOF],
    g_rel_k0: [Matrix4<f64>; N_FRAMES],
    h_inv: MMat,
    prev_iters: usize,
    pub t: f64,
}

impl SimState {
    pub fn new(cfg: SolverConfig, gains: Gains, q0: [f64; N_DOF], qdot0: [f64; N_DOF]) -> Self {
        // p_0 = M(q0)*qDot0 is only ever used through this zero case in the port (see
        // reference/math_spec.md Sec.4); the example always starts from rest.
        assert!(qdot0.iter().all(|v| *v == 0.0), "nonzero initial velocity not supported");

        let g_rel_0 = kinematics::all_g_rel(&q0);
        let g_0 = kinematics::fk(&g_rel_0);
        let j_0 = kinematics::jacobian(&q0, &g_rel_0);
        let grav_0 = gravity_wrench(&g_0);

        let a0 = 0.5;
        let mut f_gen_0 = QVec::zeros();
        for i in 0..N_DOF {
            f_gen_0[i] = cfg.h * (1.0 - a0) * gains.c_sys[i] * (q0[i] - gains.q_ref[i]);
        }
        let mut f_frame_0 = [Vec6::zeros(); N_FRAMES];
        for i in 0..N_FRAMES {
            f_frame_0[i] = grav_0[i] * (cfg.h * (1.0 - a0));
        }
        let frozen_0 = StepFrozen {
            g_rel_k: g_rel_0,
            j: j_0,
            f_gen: f_gen_0,
            f_frame: f_frame_0,
            damping_coeff: 1.0 - a0,
            q_k: q0,
        };
        let mut h_inv = MMat::identity();
        let (q1, _iters) = solve_del(&frozen_0, q0, &gains, &mut h_inv, true, true, &cfg);

        Self {
            cfg,
            gains,
            q_k0: q0,
            q_k: q1,
            g_rel_k0: g_rel_0,
            h_inv,
            prev_iters: usize::MAX, // force a fresh Jacobian on the first general step
            t: cfg.h,
        }
    }

    /// Change PD gains / target live; forces a Jacobian-approximation refresh on the next step
    /// (cheap, and avoids ever iterating against a stale `dSys`-dependent preconditioner).
    pub fn set_gains(&mut self, gains: Gains) {
        self.gains = gains;
        self.prev_iters = usize::MAX;
    }

    pub fn q(&self) -> [f64; N_DOF] {
        self.q_k
    }

    pub fn q_dot(&self) -> [f64; N_DOF] {
        core::array::from_fn(|i| (self.q_k[i] - self.q_k0[i]) / self.cfg.h)
    }

    /// Advance exactly one `h`-sized variational-integrator step.
    pub fn step(&mut self) {
        let (frozen, g_rel_k) = build_frozen_general(&self.q_k0, &self.q_k, &self.g_rel_k0, &self.gains, self.cfg.h);

        let mut guess = [0.0; N_DOF];
        for i in 0..N_DOF {
            guess[i] = 2.0 * self.q_k[i] - self.q_k0[i];
        }

        let must_recompute = self.prev_iters > self.cfg.jac_iter_threshold;
        let (q_k1, iters) = solve_del(&frozen, guess, &self.gains, &mut self.h_inv, must_recompute, false, &self.cfg);

        self.prev_iters = iters;
        self.g_rel_k0 = g_rel_k;
        self.q_k0 = self.q_k;
        self.q_k = q_k1;
        self.t += self.cfg.h;
    }

    /// Absolute world pose of every frame at the current configuration (for rendering).
    pub fn frame_poses(&self) -> [Matrix4<f64>; N_FRAMES] {
        kinematics::fk(&kinematics::all_g_rel(&self.q_k))
    }
}

pub struct Trajectory {
    pub t: Vec<f64>,
    pub q: Vec<[f64; N_DOF]>,
    pub q_dot: Vec<[f64; N_DOF]>,
}

/// Batch-run the full example (used by the native validation test): identical physics to
/// `SimState`, just driven to completion and recorded at every step.
pub fn simulate(t_end: f64, q0: [f64; N_DOF], qdot0: [f64; N_DOF], cfg: &SolverConfig) -> Trajectory {
    let n_steps = (t_end / cfg.h).round() as usize;
    let mut traj = Trajectory {
        t: Vec::with_capacity(n_steps + 1),
        q: Vec::with_capacity(n_steps + 1),
        q_dot: Vec::with_capacity(n_steps + 1),
    };
    traj.t.push(0.0);
    traj.q.push(q0);
    traj.q_dot.push(qdot0);

    let mut sim = SimState::new(*cfg, Gains::example_default(), q0, qdot0);
    traj.t.push(sim.t);
    traj.q.push(sim.q());
    traj.q_dot.push(sim.q_dot());

    for _ in 1..n_steps {
        sim.step();
        traj.t.push(sim.t);
        traj.q.push(sim.q());
        traj.q_dot.push(sim.q_dot());
    }
    traj
}
