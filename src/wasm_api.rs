//! Thin wasm-bindgen boundary: JS drives the sim with plain `Float64Array`s, no serde/JSON.

use crate::dynamics::{Gains, SimState, SolverConfig};
use crate::system_data::N_DOF;
use nalgebra::UnitQuaternion;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Sim {
    state: SimState,
}

fn gains_from_js(p: &[f64], d: &[f64], q_des_rad: &[f64]) -> Gains {
    let mut p4 = [0.0; 4];
    let mut d4 = [0.0; 4];
    let mut q4 = [0.0; 4];
    p4.copy_from_slice(&p[0..4]);
    d4.copy_from_slice(&d[0..4]);
    q4.copy_from_slice(&q_des_rad[0..4]);
    Gains::new(p4, d4, q4)
}

#[wasm_bindgen]
impl Sim {
    /// `p_gain`, `d_gain`, `q_des_rad`: length-4 arrays for the 4 PD-controlled joints
    /// (base yaw, shoulder pitch, elbow pitch, instrument roll). Starts from rest at q=0,
    /// matching the MATLAB example.
    #[wasm_bindgen(constructor)]
    pub fn new(p_gain: &[f64], d_gain: &[f64], q_des_rad: &[f64]) -> Sim {
        let gains = gains_from_js(p_gain, d_gain, q_des_rad);
        let state = SimState::new(SolverConfig::default(), gains, [0.0; N_DOF], [0.0; N_DOF]);
        Sim { state }
    }

    pub fn set_gains(&mut self, p_gain: &[f64], d_gain: &[f64], q_des_rad: &[f64]) {
        self.state.set_gains(gains_from_js(p_gain, d_gain, q_des_rad));
    }

    /// Fixed physics timestep (seconds); drive `step()` at this rate (e.g. an accumulator in
    /// the render loop) for results matching the reference trajectory.
    pub fn h(&self) -> f64 {
        self.state.cfg.h
    }

    pub fn t(&self) -> f64 {
        self.state.t
    }

    pub fn step(&mut self) {
        self.state.step();
    }

    pub fn q(&self) -> Vec<f64> {
        self.state.q().to_vec()
    }

    /// Flat `[x,y,z, qx,qy,qz,qw] * N_FRAMES` world-frame pose array, ready for
    /// `Object3D.position.set` / `Object3D.quaternion.set` in Three.js.
    pub fn frame_poses(&self) -> Vec<f64> {
        let poses = self.state.frame_poses();
        let mut out = Vec::with_capacity(poses.len() * 7);
        for g in poses.iter() {
            let r = g.fixed_view::<3, 3>(0, 0).into_owned();
            let p = g.fixed_view::<3, 1>(0, 3);
            let q = UnitQuaternion::from_matrix(&r);
            let c = q.quaternion().coords; // nalgebra order: (x, y, z, w), matches Three.js
            out.push(p[0]);
            out.push(p[1]);
            out.push(p[2]);
            out.push(c[0]);
            out.push(c[1]);
            out.push(c[2]);
            out.push(c[3]);
        }
        out
    }
}
