use crate::se3::{self, Vec6};
use crate::system_data::*;
use nalgebra::{Matrix4, SMatrix, Vector3};

pub type Jac = SMatrix<f64, 6, N_DOF>;

fn mat4_from_flat_row_major(f: &[f64; 16]) -> Matrix4<f64> {
    // f is row-major (matches MATLAB g_ref(i,j) order emitted by gen_rust_data.py)
    Matrix4::new(
        f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8], f[9], f[10], f[11], f[12], f[13],
        f[14], f[15],
    )
}

/// Relative transform g_rel(i) for frame `i` given the full 16-vector of generalized coords `q`.
pub fn g_rel(i: usize, q: &[f64; N_DOF]) -> Matrix4<f64> {
    let g_ref = mat4_from_flat_row_major(&G_REF[i]);
    match JOINT_TYPE[i] {
        1 => {
            let qi = q[Q_INDICES[i][0] as usize];
            let axis = Vec6::from_row_slice(&SCREW_AXIS[i]);
            g_ref * se3::exp_se3_screw(&axis, qi)
        }
        2 => {
            let a = Q_INDICES[i][0] as usize;
            let b = Q_INDICES[i][1] as usize;
            let l = SEG_LEN[i];
            let mut xi = Vec6::from_row_slice(&XI_C[i]);
            xi[0] += q[a]; // Ba column 0 = e1 (bending about local x)
            xi[1] += q[b]; // Ba column 1 = e2 (bending about local y)
            g_ref * se3::cay_se3(&(xi * l))
        }
        _ => unreachable!(),
    }
}

/// Forward kinematics: absolute world pose of every frame, for a given relative-transform set.
pub fn fk(g_rel_all: &[Matrix4<f64>; N_FRAMES]) -> [Matrix4<f64>; N_FRAMES] {
    let mut g = [Matrix4::identity(); N_FRAMES];
    for i in 0..N_FRAMES {
        let parent_g = if PARENT[i] < 0 {
            Matrix4::identity()
        } else {
            g[PARENT[i] as usize]
        };
        g[i] = parent_g * g_rel_all[i];
    }
    g
}

pub fn all_g_rel(q: &[f64; N_DOF]) -> [Matrix4<f64>; N_FRAMES] {
    core::array::from_fn(|i| g_rel(i, q))
}

/// Body-fixed Jacobian (6 x N_DOF) of every frame, evaluated at `q` (own-column strain) with
/// the already-computed relative transforms `g_rel_all` (own-column AdInv propagation).
pub fn jacobian(q: &[f64; N_DOF], g_rel_all: &[Matrix4<f64>; N_FRAMES]) -> [Jac; N_FRAMES] {
    let mut j: [Jac; N_FRAMES] = [Jac::zeros(); N_FRAMES];
    for i in 0..N_FRAMES {
        let mut ji = if PARENT[i] < 0 {
            Jac::zeros()
        } else {
            se3::ad_inv(&g_rel_all[i]) * j[PARENT[i] as usize]
        };
        match JOINT_TYPE[i] {
            1 => {
                let qidx = Q_INDICES[i][0] as usize;
                let axis = Vec6::from_row_slice(&SCREW_AXIS[i]);
                ji.set_column(qidx, &axis);
            }
            2 => {
                let a = Q_INDICES[i][0] as usize;
                let b = Q_INDICES[i][1] as usize;
                let l = SEG_LEN[i];
                let mut xi = Vec6::from_row_slice(&XI_C[i]);
                xi[0] += q[a];
                xi[1] += q[b];
                let t = se3::cay_rtd_se3(&(-(xi * l)));
                let col_a = l * t.column(0);
                let col_b = l * t.column(1);
                ji.set_column(a, &col_a);
                ji.set_column(b, &col_b);
            }
            _ => unreachable!(),
        }
        j[i] = ji;
    }
    j
}

/// World "up" unit vector used by the gravity term, body-frame projected per frame via R(i)'.
pub fn gravity_dir() -> Vector3<f64> {
    Vector3::new(0.0, 0.0, 1.0)
}
