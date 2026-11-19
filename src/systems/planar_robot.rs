//! `systemDef_planarNLinkPendulum("nLinks", 2, "d", 0)` - exact numbers probed live from the
//! installed MATLAB toolbox (reference/ocp_planar_frame_data.json), not re-derived by hand.

use crate::rigid_chain::ChainDef;
use crate::se3::vec6;
use nalgebra::{Matrix4, Vector3};

pub const N: usize = 2;

fn se3_identity_x(x: f64) -> Matrix4<f64> {
    let mut g = Matrix4::<f64>::identity();
    g[(0, 3)] = x;
    g
}

pub fn def() -> ChainDef<N> {
    let axis = vec6(Vector3::new(0.0, 1.0, 0.0), Vector3::new(0.0, 0.0, -0.5));
    ChainDef {
        screw_axis: [axis, axis],
        g_ref: [se3_identity_x(0.5), se3_identity_x(1.0)],
        parent: [-1, 0],
        m_gen_diag: [
            [5.8946275219220495e-5, 0.083362806470942935, 0.083362806470942935, 1.0, 1.0, 1.0],
            [5.8946275219220495e-5, 0.083362806470942935, 0.083362806470942935, 1.0, 1.0, 1.0],
        ],
        mass: [1.0, 1.0],
        g_const: 9.81,
    }
}
