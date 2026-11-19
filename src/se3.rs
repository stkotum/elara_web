//! SE(3)/SO(3) Lie group operations, ported 1:1 from ELARA's `elara/math/*.m`.
//! Twist convention: xi = [omega(3); v(3)] (angular first), matching ELARA's `skew.m`.

use nalgebra::{Matrix3, Matrix4, Matrix6, SMatrix, Vector3};

pub type Vec6 = SMatrix<f64, 6, 1>;

pub fn skew_so3(w: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -w.z, w.y, w.z, 0.0, -w.x, -w.y, w.x, 0.0)
}

pub fn omega_v(xi: &Vec6) -> (Vector3<f64>, Vector3<f64>) {
    (Vector3::new(xi[0], xi[1], xi[2]), Vector3::new(xi[3], xi[4], xi[5]))
}

pub fn vec6(omega: Vector3<f64>, v: Vector3<f64>) -> Vec6 {
    Vec6::new(omega.x, omega.y, omega.z, v.x, v.y, v.z)
}

pub fn se3(r: Matrix3<f64>, p: Vector3<f64>) -> Matrix4<f64> {
    let mut g = Matrix4::identity();
    g.fixed_view_mut::<3, 3>(0, 0).copy_from(&r);
    g.fixed_view_mut::<3, 1>(0, 3).copy_from(&p);
    g
}

pub fn rot(g: &Matrix4<f64>) -> Matrix3<f64> {
    g.fixed_view::<3, 3>(0, 0).into()
}

pub fn pos(g: &Matrix4<f64>) -> Vector3<f64> {
    g.fixed_view::<3, 1>(0, 3).into()
}

/// SE3Matrix.m inverse: gInv = [R'; -R'x; 0 0 0 1]
pub fn se3_inv(g: &Matrix4<f64>) -> Matrix4<f64> {
    let r = rot(g);
    let p = pos(g);
    let rt = r.transpose();
    se3(rt, -(rt * p))
}

/// Adjoint of g (6x6), `SE3.Ad`: [[R,0],[skew(x)R, R]]
pub fn ad(g: &Matrix4<f64>) -> Matrix6<f64> {
    let r = rot(g);
    let p = pos(g);
    let mut m = Matrix6::zeros();
    m.fixed_view_mut::<3, 3>(0, 0).copy_from(&r);
    m.fixed_view_mut::<3, 3>(3, 0).copy_from(&(skew_so3(&p) * r));
    m.fixed_view_mut::<3, 3>(3, 3).copy_from(&r);
    m
}

/// `lAdSE3Inv.m`: Ad(g^-1) in closed form: [[R',0],[-R'skew(x), R']]
pub fn ad_inv(g: &Matrix4<f64>) -> Matrix6<f64> {
    let r = rot(g);
    let p = pos(g);
    let rt = r.transpose();
    let mut m = Matrix6::zeros();
    m.fixed_view_mut::<3, 3>(0, 0).copy_from(&rt);
    m.fixed_view_mut::<3, 3>(3, 0).copy_from(&(-(rt * skew_so3(&p))));
    m.fixed_view_mut::<3, 3>(3, 3).copy_from(&rt);
    m
}

/// `caySO3.m`: R = I + 4/(4+om'om) * (omHat + 1/2 omHat^2)
pub fn cay_so3(om: &Vector3<f64>) -> Matrix3<f64> {
    let k = om.dot(om);
    let oh = skew_so3(om);
    let i3 = Matrix3::<f64>::identity();
    i3 + (4.0 / (4.0 + k)) * (oh + 0.5 * oh * oh)
}

/// `caySE3.m` ([Dem+14] eq.19)
pub fn cay_se3(xi: &Vec6) -> Matrix4<f64> {
    let (om, v) = omega_v(xi);
    let k = om.dot(&om);
    let oh = skew_so3(&om);
    let i3 = Matrix3::<f64>::identity();
    let r = cay_so3(&om);
    let x = (4.0 / (4.0 + k)) * (i3 + 0.5 * oh + 0.25 * om * om.transpose()) * v;
    se3(r, x)
}

/// `cayInvSE3.m` / `SE3.cayInv`
pub fn cay_inv_se3(g: &Matrix4<f64>) -> Vec6 {
    let r = rot(g);
    let p = pos(g);
    let i3 = Matrix3::<f64>::identity();
    let omega_hat = (2.0 / (1.0 + r.trace())) * (r - r.transpose());
    let omega = Vector3::new(omega_hat[(2, 1)], omega_hat[(0, 2)], omega_hat[(1, 0)]);
    let inv = (r + i3).try_inverse().expect("R+I singular in cayInvSE3");
    let v = 2.0 * inv * p;
    vec6(omega, v)
}

/// `cayRTDSE3.m` (T operator, right-trivialized derivative of cay)
pub fn cay_rtd_se3(xi: &Vec6) -> Matrix6<f64> {
    let (om, v) = omega_v(xi);
    let k = om.dot(&om);
    let oh = skew_so3(&om);
    let vh = skew_so3(&v);
    let i3 = Matrix3::<f64>::identity();
    let mut t = Matrix6::zeros();
    t.fixed_view_mut::<3, 3>(0, 0)
        .copy_from(&((2.0 / (4.0 + k)) * (2.0 * i3 + oh)));
    t.fixed_view_mut::<3, 3>(3, 0)
        .copy_from(&((1.0 / (4.0 + k)) * vh * (2.0 * i3 + oh)));
    t.fixed_view_mut::<3, 3>(3, 3)
        .copy_from(&(i3 + (1.0 / (4.0 + k)) * (2.0 * oh + oh * oh)));
    t
}

/// `cayRTDInvSE3.m` (Tinv operator, inverse right-trivialized derivative of cay)
pub fn cay_rtd_inv_se3(xi: &Vec6) -> Matrix6<f64> {
    let (om, v) = omega_v(xi);
    let oh = skew_so3(&om);
    let vh = skew_so3(&v);
    let i3 = Matrix3::<f64>::identity();
    let a = i3 - 0.5 * oh; // reused in two blocks
    let mut tinv = Matrix6::zeros();
    tinv.fixed_view_mut::<3, 3>(0, 0)
        .copy_from(&(i3 - 0.5 * oh + 0.25 * om * om.transpose()));
    tinv.fixed_view_mut::<3, 3>(3, 0).copy_from(&(-0.5 * a * vh));
    tinv.fixed_view_mut::<3, 3>(3, 3).copy_from(&a);
    tinv
}

/// `expSE3Screw.m`: constant-screw-axis exponential, [MLS94] eq.2.36.
/// xi = [omega; v], |omega| assumed 1 (or omega==0 for pure translation).
pub fn exp_se3_screw(xi: &Vec6, theta: f64) -> Matrix4<f64> {
    let (om, v) = omega_v(xi);
    let i3 = Matrix3::<f64>::identity();
    if om.norm() < 1e-12 {
        return se3(i3, v * theta);
    }
    let oh = skew_so3(&om);
    let r = i3 + theta.sin() * oh + (1.0 - theta.cos()) * (oh * oh);
    let p = (i3 - r) * (om.cross(&v)) + om * (om.dot(&v) * theta);
    se3(r, p)
}
