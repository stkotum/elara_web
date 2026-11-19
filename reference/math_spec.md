# ELARA "simulation_rigid_flexible_robot_PD_control" - Implementable Math Spec

Source: installed MATLAB Add-On at
`C:\Users\go79ran\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\ELARA Robot Simulation and Control Toolbox\`
(NOT the GitHub main branch, which uses an older `+elara` namespace with different names).
Verified by actually running the example via `matlab -batch` (nDoF=16, solver converges,
residuals ~1e-13). All formulas below are traceable to source files cited inline.

## 1. SE(3)/SO(3) conventions

- Twist vector order everywhere: `xi = [omega; v]` (angular first, 3), (linear second, 3).
  Source: `elara/math/skew.m`, `SE3.m`. Angular/linear, NOT linear/angular.
- `se(3)` 4x4 matrix (hat map, `skew.m`/`skewSE3.m`):
  `skew(xi) = [ skewSO3(omega), v ; 0 0 0, 0 ]`, with `skewSO3([a,b,c]) = [0 -c b; c 0 -a; -b a 0]`.
- `SE3` group element stored as 4x4 homogeneous matrix `g = [R, x; 0 0 0, 1]` (`SE3Matrix.m`).
  Inverse: `gInv = [R'; -R'x ; 0 0 0 1]` (`invSE3Matrix.m`).
- Adjoint `Ad(g)` (6x6, matches `SE3.Ad` method in `SE3.m`):
  ```
  Ad(g) = [ R,          0 ]
          [ skew(x)*R,  R ]
  ```
  `AdInv(g) = Ad(g^-1)` (closed form in `lAdSE3Inv.m`, matches `[MLS94, p.56]`):
  ```
  AdInv(g) = [ R',            0  ]
             [ -R'*skew(x),   R' ]
  ```
  Code calls this operator `lAdSE3Inv` ("large Ad, inverse").
- **Two retraction maps exist in the library: exponential (`expSE3.m`) and Cayley (`caySE3.m`).
  The variational integrator used by this example (`MBSimIntegratorVarIntBroyden` /
  `integrateMBSDynamics_Broyden.m`) uses the CAYLEY map exclusively**, not the exponential map,
  for both the flexible-joint kinematics and the discrete velocity reconstruction. The
  exponential map (`expSE3Screw.m`, a constant-screw-axis variant from [MLS94]) is used only
  for the four 1-DOF *rigid revolute joints* (jointType==1), which is a closed-form pure
  rotation since their screw axis has zero pitch (see Sec. 4).
- Cayley map for SE(3) (`caySE3.m`, source [Dem+14] eq.19), closed form:
  ```
  om = xi(1:3); v = xi(4:6); omH = skew(om)
  R = caySO3(om)                                          # caySO3.m: R = I + 4/(4+om'om) * (omH + 1/2 omH^2)
  x = 4/(4 + om'om) * ( I + 1/2 omH + 1/4 om*om' ) * v
  ```
- Inverse Cayley map for SE(3) (`cayInvSE3.m` / `SE3.cayInv`), source [Dem+14] p.10:
  ```
  omegaHat = 2/(1+trace(R)) * (R - R')         # vee(omegaHat) = omega
  v = 2 * (R + I)^-1 * p                        # p = translation part of g
  xi = [omega; v]
  ```
- Right-trivialized derivative of Cayley map, `dcay`/`cayRTDSE3.m` (source [Kob14] eq.17,
  [Dem+14] p.11), 6x6 linear operator T(xi) such that it maps a velocity in se(3) through the
  Cayley retraction's tangent map:
  ```
  om = xi(1:3); v = xi(4:6); omH = skew(om); vH = skew(v); k = om'*om
  T[1:3,1:3] = 2/(4+k) * (2*I + omH)
  T[1:3,4:6] = 0
  T[4:6,1:3] = 1/(4+k) * vH * (2*I + omH)
  T[4:6,4:6] = I + 1/(4+k) * (2*omH + omH^2)
  ```
- Inverse right-trivialized derivative, `cayRTDInvSE3.m` (source [KM11] eq.34, [Dem+14] eq.20):
  ```
  om = xi(1:3); v = xi(4:6); omH = skew(om); vH = skew(v)
  Tinv[1:3,1:3] = I - 1/2*omH + 1/4*(om*om')
  Tinv[1:3,4:6] = 0
  Tinv[4:6,1:3] = -1/2 * (I - 1/2*omH) * vH
  Tinv[4:6,4:6] = I - 1/2*omH
  ```
  This is the one actually used at every integration step (inertia term + Broyden Jacobian),
  so it is the single most important 6x6 operator to port correctly.

## 2. The robot (exact numbers, `systemDef_rigid_flexible_robot.m` + `systemDef_rigid_robot.m`)

Serial open chain, 4 links, base-to-tip order. All links parented to the previous one
(parentLink = 0,1,2,3). `g0 = eye(4)` (world = base frame, confirmed by probe: `g0` default).

**Link 1 (rigid, revolute, base):**
- Joint axis Z = `[0,0,1,0,0,0]` (rotation about local z, zero pitch, so pure rotation, `v=0`)
- `l1 = 0.25 m`, cylinder radius `r=0.05 m`, density `rho=2700 kg/m^3` (aluminium)
- `m1 = r^2*pi*rho*l1 = 5.3014 kg` (matches probe output exactly)
- `I_sym = 1/2*m*r^2` (about cylinder's own long axis), `I_long = 1/4*m*r^2 + 1/12*m*l^2`
  (about transverse axes), giving `J1 = diag([I_long, I_long, I_sym])`
- `g_J1_COM = SE3(I, [0,0,l1/2])` (COM at half the link length along local z)
- Joint damping `d = 0.01`, no explicit joint stiffness here (stiffness set later via `cSys`)

**Link 2 (rigid, revolute):**
- Constant joint rotation `g_J_const = SE3([0 -1 0; 0 0 -1; 1 0 0], 0)` applied between link1's
  tip frame and link2's joint frame (rotates the joint axis 90 deg, i.e. elbow axis is
  perpendicular to the base yaw axis, a standard "yaw then pitch" robot arm layout)
- `l2 = 0.4 m`, radius `0.03 m`, same rho=2700, giving `m2 = 3.0536 kg` (matches probe)
- `J2 = diag([I_sym, I_long, I_long])` (note axis order differs from link 1: link extends
  along local x here, not z: `g_J1_J2 = SE3(I,[l,0,0])`, `g_J1_COM = SE3(I,[l/2,0,0])`)
- `d = 0.01`

**Link 3 (rigid, revolute):** identical formulas to link 2 with `l3 = 0.4 m`, giving
  `m3 = 3.0536 kg` (matches probe). `g_J_const = eye(4)` (no extra rotation vs link 2's frame).

**Link 4 (rigid 1-DOF joint + flexible beam, the "surgical instrument"):**
- `parentLink = 3`. First frame of this link is itself a screw joint (NOT part of the beam):
  `jointAxis = [1,0,0,0,0,0]` (rotation about local x). `g_J_B = SE3(I,[0,0,0.1])`,
  `g_ref = SE3([0 0 1; 0 1 0; -1 0 0], [l3/2+0.05, 0, 0])` relative to link 3's frame.
  `d = opts.dJoints = 1e-2` (same joint damping as the rigid joints).
- This 4th screw joint, together with the 3 rigid links above, is exactly the
  `nLinksRigid = 4` referenced in the example script (`qDes` has 4 entries, one per screw
  joint; `cSys`/`dSys` for these 4 DOF get overwritten with the PD gains, see Sec. 5).
- Beam: `nSeg = 6` segments (script does not override the default), `L = 0.4 m` total length,
  giving segment length `l = L/nSeg = 0.4/6 = 0.066667 m` (matches probe frame table exactly).
  Beam is NOT a cantilever (`isCantilever=false`): since its own first frame is the screw
  joint above, the beam itself occupies `nSeg+1 = 7` node frames (global frames 5..10 from
  the probe; segment-joint, i.e. jointType==2, count = 6, matching `2*nSeg = 12` DOF below).
- Allowed deformation modes `Ba = [eye(2); zeros(4,2)]`: **only the first two strain
  components (bending about local x and y) are free generalized coordinates per segment.**
  Constrained modes `Bc = [zeros(2,4); eye(4)]` fix torsion (omega_z), both shear
  components (v_x, v_y), and axial stretch (v_z) to the reference strain
  `xiRef = [0;0;0;0;0;1]` per segment (i.e. straight, unit-stretch, untwisted, unsheared
  reference configuration, a pure-bending Euler-Bernoulli-like reduced beam model, NOT a
  full 6-DOF-per-segment geometrically exact rod). This gives 2 DOF/segment * 6 segments = 12 DOF.
- Total `nDoF = 4 (screw joints) + 12 (beam bending DOF) = 16`. **Confirmed empirically**
  (`MBSim.MBSys.nDoF == 16` when run in MATLAB).
- Beam material: `beamParams_ASA_round(radius=2.75e-3)`, circular cross-section,
  `radius = 2.75 mm`, `A = pi*r^2`, `I_x = I_y = pi/4*r^4`, `J_P (polar) = pi/2*r^4`.
  Material = "ASA" 3D-printed plastic: `rho = 1070 kg/m^3`, `E = 2279.5e6 Pa` (default, NOT
  overridden by the example), `nu = 0.3459`, giving `G = E/(2*(1+nu))`.
  Damping override: `linkF.beamPars.d = opts.dBeam = ones(6,1)*1e-3` (example's actual call
  passes `"dBeam", ones(6,1)*1e-3`, i.e. uniform Kelvin-Voigt-type linear strain-rate damping
  of 1e-3 on all 6 strain components, though only the 2 "allowed" (Ba) components actually end
  up in `dSys` (see Sec. 3/5 formula `dSys = l * Ba' * beamPars.d`).
- Two point masses glued to the beam (steel balls, rho=8211 kg/m^3):
  - At the first beam node (the joint-side end): `rBall_1 = 0.03 m`, giving
    `m_b_1 = rho*4/3*pi*rBall_1^3 = 0.92864 kg` (probe: frame4 `m_a=0.92864`),
    `J_b_1 = I * m_b_1 * 2/5 * rBall_1^2` (solid sphere inertia).
  - At the last beam node (tip/end-effector): `rBall_2 = 0.007 m`, giving
    `m_b_2 = 0.011797 kg` (probe: frame10 `m_a=0.011797`), `J_b_2` likewise.
  - These add directly into that frame's `MGen` (6x6) and `m` on top of the beam segment's own
    half-segment mass contribution (see Sec. 3 mass assembly: `factors=[0.5,1,1,1,1,1,0.5]`
    for the 7 beam nodes, i.e. trapezoidal lumping along the segment chain).
  - TCP (tool center point) = last beam node, `g_B_TCP = eye(4)` offset.
- Gravity: `simPars.g = 9.81` (default, scalar magnitude), direction/sign: see Sec. 6 (flagged
  as needing empirical sign validation).

## 3. Mass matrix and forward kinematics assembly (`assembleMBSystem.m`, `MBSystemNum.m`)

**Frame graph.** Each rigid link = 1 frame. Each flexible (non-cantilever) link = `nSeg+1`
frames: frame 1 of the link = the screw joint connecting it to its parent (jointType=1),
frames 2..nSeg+1 = the beam segment nodes (jointType=2). Frames are numbered globally,
base to tip, in link order (confirmed by probe: frames 1-4 = the four screw joints, frames
5-10 = the six beam-segment joints).

**Per-frame generalized mass `MGen` (6x6, block-diag [J_rot, m*I_3] in general, but may be a
full symmetric 6x6 if off-diagonal coupling exists):**
- Rigid link frame: `MGen = blkdiag(link.J, link.m*eye(3))`.
- Flexible beam node frame: `MGen = l * factor(node) * beamPars.Mgen (+ external point mass
  M_a(node) if present)`, where `beamPars.Mgen = blkdiag(beamPars.J, beamPars.m*eye(3))`,
  `beamPars.J = rho*diag([I_x, I_y, J_P])`, `beamPars.m = rho*A` (mass per unit length), and
  `factor = [0.5, 1,1,1,1,1, 0.5]` for the `nSeg+1=7` nodes of a non-cantilever beam
  (trapezoidal/lumped-mass rule: end nodes get half a segment's worth of mass, interior
  nodes get a full segment's worth; note both neighboring segments would contribute to an
  interior node in a strict trapezoidal rule, but this code's `factor` array as written gives
  interior nodes weight 1 x l, not (l_left+l_right)/2 x ...; with uniform l these are the same
  thing, so no discrepancy for this example since all segments share the same length).
- `m(frame) = l*factor*beamPars.m (+ m_a)`.

**Joint transformation `g_rel` per frame (`computeJointTransformations`):**
- jointType 1 (screw): `g_rel = g_ref * expSE3Screw(X, q_i)` where `X = AdInv(g_J_B) * jointAxis`
  (joint axis is transformed from the joint-definition frame into the body/CoM frame via the
  inverse adjoint of `g_J_B`), `q_i` = the single generalized coordinate (joint angle) for that
  frame.
- jointType 2 (flexible): `xi = Ba(frame) * q_i + xiC(frame)`, `xiC = Bc*Bc' * xiRef(segment)`
  (projects the reference strain onto the constrained subspace), then
  `g_rel = caySE3( xi * l )` (segment length `l` scales the reduced strain `xi` into the
  actual SE(3) relative transform over that segment length, the Cayley map of the *scaled*
  strain vector).

**Forward kinematics:** `g(1) = g0 * g_rel(1)`; `g(i) = g(parent(i)) * g_rel(i)` for i>1,
recursively frame-by-frame from base to tip (`computeFwdKinFast`).

**Geometric Jacobian `J` (6 x nDoF x nFrames), built incrementally (`computeGeomJacobianFast`):**
- Own-column block of frame `i` (`ii==iFrm`):
  - jointType 1: `J(:,qIdx,i) = X(i)` (the screw axis itself, constant per frame).
  - jointType 2: `J(:,qIdx,i) = l(i) * cayRTDSE3(-xi(i)*l(i)) * Ba(i)` (right-trivialized
    Cayley tangent map of the NEGATIVE scaled strain, times the segment length, times the
    mode-selection matrix).
- Ancestor-column blocks (`ii<iFrm`, `ii` an ancestor of `iFrm`):
  `J(:,qIdx_ii,iFrm) = AdInv(g_rel(iFrm)) * J(:,qIdx_ii,parent(iFrm))`, i.e. propagate the
  parent frame's Jacobian columns through the inverse adjoint of the current relative
  transform (standard body-fixed/right-trivialized recursive Jacobian for a serial chain).

**Mass matrix:** `M(q) = sum_i J(:,:,i)' * MGen(i) * J(:,:,i)` over all frames (classic
Jacobian-weighted sum of per-body spatial inertias, expressed in generalized coordinates).

**Discrete (interval) body-fixed velocity `eta_k` per frame, for the step from `k` to `k+1`**
(`computeDiscreteAbsoluteVelocities`), computed recursively base to tip:
```
eta_k(1)   = cayInvSE3( invSE3Matrix(g_rel_k(1)) * g_rel_k1(1) ) / h
eta_k(i>1) = cayInvSE3( invSE3Matrix(g_rel_k(i)) * caySE3(h*eta_k(i-1)) * g_rel_k1(i) ) / h
```
Note the recursion literally uses `eta_k(i-1)`, i.e. it assumes frames are indexed in strict
serial chain order matching physical parent-child adjacency (true for this example's serial
topology; would need generalizing to `eta_k(parent(i))` for a branching tree, not needed here).

## 4. Variational integrator: DEL residual and Broyden solve

Source: `integrateMBSDynamics_Broyden.m` + `solveImplicitDELEquBroyden.m`. This is a Lie-group
discrete Euler-Lagrange (DEL) scheme on SE(3)^nFrames (per [Lee08]-style Lie group variational
integrators), solved implicitly with Broyden's "good" quasi-Newton method instead of full
Newton (because the exact Jacobian of the DEL residual is expensive; Broyden amortizes it).

**Per-step unknown:** `q_k1` (the coordinates at the *next* time step; `q_k`, `q_k0`, current
and previous, are already known). One step advances `(q_k0, q_k) -> q_k1`.

**Residual function** (`MBSDynamics_DEL_implicitFun`, local function in
`solveImplicitDELEquBroyden.m`), given `a` = the generalized-trapezoidal weight for THIS call
(`a=1/2` fixed for the very first step only; `a = aTrapez` solver setting, which the example
sets to `0`, for every subsequent step):
```
q_k_diff = q_k1 - q_k                      # forward difference, the unknown's own step
DEL_res = f_gen_k + (1-a) * dSys .* q_k_diff            # (if any dSys != 0)
for each frame i:
    DEL_res += JGeom(:,:,i)' * ( cayRTDInvSE3(eta_k(i)*h)' * MGen(i) * eta_k(i)  +  f_frame_k(:,i) )
```
`JGeom` is evaluated ONCE per step at `q_k` (not re-evaluated each Broyden iteration, a
"frozen Jacobian" approximation that's part of why Broyden, not full Newton, is adequate).
`eta_k` is recomputed each iteration from the current `q_k1` guess via
`computeDiscreteAbsoluteVelocities(g_rel_k, g_rel_k1(q_k1), h)`.

**`f_gen_k` (generalized force input to the residual), assembled once per step in
`integrateMBSDynamics_Broyden.m` (this is where joint stiffness/PD and external inputs enter):**
```
f_gen_k = h * cSys .* (q_k - qRef)  -  h * computeInputMatrixFast(g_rel_k) * u(:,k)  +  aLoop * dSys .* (q_k - q_k0)
```
(`aLoop = solverConfig.aTrapez`; this term is EXTRA to the `(1-a)*dSys*q_k_diff` term already
inside the residual function above, with `aTrapez=0` as in the example, this extra term
vanishes, but the `(1-a)*dSys*(q_k1-q_k)` term inside the residual does NOT vanish since
`a` there equals `aLoop=0` too, giving `(1-0)=1`, i.e. full one-sided/backward-difference
linear damping `+dSys.*(q_k1-q_k)` is active). **Net effect for this example (aTrapez=0): the
PD damping term acting on the system is `dSys .* (q_k1 - q_k)`, a first-order/one-sided
discrete derivative, NOT a centered/trapezoidal one.** This is a subtlety worth a unit test
when porting (verify against the MATLAB-generated q(t) that damping is present and matches
this one-sided form, not zero and not the centered form).

**`f_frame_k` (per-frame body-fixed wrench: inertia, gravity, external), assembled per step:**
```
f_frame_k_b = -h*f_frame_k_b_external + h*computeBodyfixedFrameForces(g_k, f_ext_spatial, MBSys, simPars) - computeInertiaTerm(MBSys, eta_k0, h)
computeInertiaTerm: f_inertia(i) = cayRTDInvSE3(-eta_k0(i)*h)' * MGen(i) * eta_k0(i)   # note NEGATIVE eta_k0*h here
```
`eta_k0` is the discrete velocity from the PREVIOUS interval (k-1,k); the inertia term uses the
PAST velocity with a *negated* argument to `cayRTDInvSE3`, contrasted with the residual's
own-interval inertia term `cayRTDInvSE3(+eta_k(i)*h)' * MGen(i) * eta_k(i)` (positive argument).
This left/right-momentum distinction (discrete Legendre transforms FL-/FL+) is the standard
structure of a Lie group variational integrator and must be reproduced exactly (sign matters).

**Gravity and external force (`computeBodyfixedFrameForces.m`):**
```
f_frame_b(i) = [ R(i)' * -f_ext_spatial_torque(i) ; R(i)' * -f_ext_spatial_force(i) ]
             + g_const * [ m_a(i) * cross(x_a(i), R(i)'*[0;0;1]) ; m(i) * R(i)'*[0;0;1] ]
```
`g_const = simPars.g = 9.81`, world gravity direction assumed to be `+z` in this formula as
written (**flag: verify sign/direction against the reference trajectory, see Sec. 6**).
`R(i)` = rotation part of the frame's ABSOLUTE configuration `g(i)` (world orientation), so
this correctly rotates the constant world-frame gravity vector into the body frame.

**Initial step (k=1 to k=2), handled specially to seed from a given initial velocity
`qDot0` via the "left" discrete Legendre transform:**
```
a = 1/2 (always, regardless of aTrapez)
p_0 = M(q0) * qDot0                                    # continuous-time initial momentum
f_frame_0_b = -h*(1-a)*f_ext + h*(1-a)*computeBodyfixedFrameForces(g_0, f_ext_s, ...)
f_gen_0     = h*(1-a)*cSys.*(q0-qRef) - h*(1-a)*computeInputMatrixFast(g_rel_0)*u(:,1) - p_0
```
then solved with the SAME `solveImplicitDELEquBroyden` (with `a=1/2` passed through), using
initial guess `q_k1_guess = q0` itself (`q_k0 = q_k = q0` for this call) and
`forceSolverIteration = true` (always iterate at least once, since residual may spuriously
read as converged at the trivial guess).

**Loop step (k=2..nSteps):** initial guess for the Broyden solve is `q_k1_guess = 2*q_k - q_k0`
(explicit extrapolation, "IG2" from Lee et al., a second-order predictor).

**Broyden's "good" method** (`solveImplicitDELEquBroyden.m`):
1. Evaluate `resDEL` at the extrapolated guess.
2. If `|resDEL| <= errorMargin` (and not forced): done, 0 iterations.
3. Else, if `updateInvJacobian` is true (true on step 1, and whenever the previous step took
   more than `JacobianIterationThreshold` iterations, default threshold in this example is
   5), recompute the *approximate* inverse Jacobian analytically as:
   ```
   MBeam = sum_i JGeom(:,:,i)' * cayRTDInvSE3(eta_k(i)*h)' * MGen(i) * JGeom(:,:,i)
   H_k = inv( MBeam/h + diag(dSys) )        # NOTE: full dSys used here (not (1-a)*dSys);
                                              # comment in source says this empirically
                                              # converges better than the "exact" (1-a) form
   ```
   (`MBeam` here is effectively the discrete mass matrix re-derived through the Cayley
   tangent map rather than the plain `computeMassMatrix`; these coincide to leading order
   but are not identical, so use this exact formula, not the continuous mass matrix, for the
   Jacobian approximation.)
4. Otherwise reuse the inverse Jacobian `H_k` carried over from the previous time step
   (this is the whole point of Broyden's method: avoid recomputing/inverting a matrix every
   iteration or even every time step).
5. Newton-like update: `q_k1 <- q_k1 - H_k * resDEL`. Re-evaluate `resDEL`.
6. If converged, stop. Else Broyden "good" rank-1 update of the inverse Jacobian:
   ```
   s_k = q_k1_new - q_k1_old         (the step taken)
   y_k = resDEL_new - resDEL_old     (the resulting residual change)
   H_k <- H_k + ( (s_k - H_k*y_k) * (s_k' * H_k) ) / (s_k' * H_k * y_k)
   ```
7. Repeat up to `maxIterations` (100 default; example doesn't override).
8. Convergence tolerance for this example: `errorMargin = 5e-12` (very tight,
   `solverConfig.errorMargin`), `JacobianIterationThreshold = 5`, `h = 2^-8 ~= 0.003906 s`.

**Loop bookkeeping:** `q_dot(k) = (q_k1-q_k)/h` is stored as an OUTPUT convenience (the finite
difference of the solved trajectory), not used internally by the solver itself beyond this.
`eta` (per-frame discrete body velocity) is the "real" internal velocity representation.

## 5. PD control term, exact mechanism

The example does NOT add an explicit continuous PD torque computation; instead it repurposes
the system's own linear-stiffness/dissipation machinery (`MBSys.cSys`, `MBSys.dSys`,
`MBSys.qRef`), which normally represents the flexible beam's elastic/damping properties, and
overwrites the first 4 entries (the 4 screw-joint DOF) with PD gains:
```matlab
nLinksRigid = 4;
MBSim.MBSys.dSys(1:nLinksRigid) = 30;     % "D" gain (damping/derivative)
MBSim.MBSys.cSys(1:nLinksRigid) = 100;    % "P" gain (stiffness/proportional)
MBSim.MBSys.qRef(1:nLinksRigid) = qDes;   % desired joint angles, deg2rad([45,-45,90,-80])
```
These then flow through exactly the same `f_gen_k = h*cSys.*(q_k-qRef) + ...` and
`(1-a)*dSys.*(q_k1-q_k)` residual terms derived in Sec. 4; i.e., for this example, the "PD
controller" IS the beam's own generalized stiffness/damping formula, just applied to the rigid
joint coordinates instead of beam strains, with gains chosen by the user rather than derived
from `beamParams`. **This means: there is no separate control law to port, just initialize
`cSys[0:4]=100`, `dSys[0:4]=30`, `qRef[0:4]=qDes` on top of the already-assembled system (which
for the beam DOF 4:16 has `cSys`/`dSys`/`qRef` set from the beam stiffness formula in Sec.3,
`cSys = l*Ba'*diag(Cgen)`, `dSys = l*Ba'*beamPars.d`, `qRef = Ba'*xiRef`), then run the same
single DEL/Broyden loop for the whole 16-DOF system.** The flexible beam entries of `dSys` use
`beamPars.d = dBeam = [1e-3]*6` (example override) projected through `Ba` (only the first 2
rows of `Ba'` survive, since `Ba=[eye(2);zeros(4,2)]`), giving small but nonzero linear strain
damping on the beam bending DOF too, separate from, and in addition to, the rigid-joint PD
damping.

`qDes = deg2rad([45, -45, 90, -80])` (4 values, one per screw joint, in the order:
base yaw, shoulder pitch, elbow pitch, instrument-roll).

## 6. Flags for empirical validation against the MATLAB reference run

These are points where either (a) a global sign/direction convention could not be pinned down
with full certainty from static reading alone, or (b) multiple numerically-equivalent-looking
formulas exist in the source and only one is exercised by this example path:

1. **Gravity sign/direction** (`computeBodyfixedFrameForces.m`): formula as written uses `+g*R'*[0;0;1]`
   with no explicit minus sign anywhere in that function or its caller. Whether this nets out
   to gravity pulling the robot toward -z or +z in world frame depends on how this "frame
   force" is combined into the residual's overall sign convention (it enters as `+f_frame_k_b`
   inside the DEL residual after being folded through `+h*computeBodyfixedFrameForces(...)`
   inside the larger `f_frame_k_b = -h*f_ext + h*computeBodyfixedFrameForces(...) - inertia`
   expression). **Action: after porting, run a trivial one-link pendulum test (or just this
   full example with all PD gains set to 0 and zero initial velocity) and visually/numerically
   confirm the robot's beam tip sags in -z under its own weight, matching a MATLAB reference
   run of the same zero-gain case. Flip the gravity vector sign in the port if it sags upward.**
2. **Rectangle-vs-trapezoidal damping term, confirmed by reading** (Sec. 4): with
   `aTrapez=0`, damping is `dSys.*(q_k1-q_k)` only (one-sided). This was traced precisely
   through the code (not guessed), but should still get a regression test: run the port with
   `cSys=0, dSys>0` on one DOF with nonzero initial velocity and confirm it decays
   monotonically at the expected rate, cross-checked against a short MATLAB reference run.
3. **MEX vs .m numerics**: `MBSimIntegratorVarIntBroyden.simulateSystem` prefers a compiled
   `integrateMBSDynamics_Broyden_mex` if present on the MATLAB path over the plain `.m` file
   we read. The `.m` source (what this spec is based on) is marked `%#codegen`, meaning the
   mex build is auto-generated FROM this exact `.m` file with no manual edits, so they should
   be numerically identical (same algorithm, just compiled), but when generating the
   reference trajectory for validation, confirm which path actually ran (the console output
   line `Starting integration (varint-broyden)...` prints regardless; checking
   `exist("integrateMBSDynamics_Broyden_mex","file")` separately, or just diffing a short run's
   output against hand-stepping the .m algorithm, removes any doubt).
4. **`MBSys.frameData.ancestors`/`parent` vs. the hard-coded `iFrm-1` recursion** in
   `computeDiscreteAbsoluteVelocities`: works correctly for this example's pure serial chain;
   a Rust port can just hard-code "previous frame in the base-to-tip ordering" rather than
   building a general tree-parent lookup, since this example has no branching.
5. Not traced in depth (not exercised by this example: inputs `u` are empty,
   `uConst`/`uSampleValues` unset so `computeInputMatrixFast(...)*u` contributes exactly zero):
   cable actuation (`computeInputMatrixFast`'s `case 2` branch, `g_cm`/cable routing). Safe to
   skip entirely for a first port of this specific example; `nInputs=4` exist (one per screw
   joint, per probe output) but are never driven (`simPars.uConst` is left empty in the
   script), so the `-h*B*u` term in `f_gen_k` is always exactly zero here.

## 7. Rendering: configuration to 3D geometry

`simRes.g` is already a `(4,4,nFrames,nSteps)` array of ABSOLUTE SE(3) poses for every frame at
every time step (world frame, computed by `computeFwdKinFast` each step), the integrator
output already IS the renderable pose data, no further forward-kinematics needed at render
time, just index by time step.

- **Rigid links** (frames 1-4): draw as a box/cylinder between the frame's own origin and
  context from `bBoxSize`/`g_bbox` in the link definition (`MBRigidLinkVisualization.m`,
  not read in full, but `bBoxSize` gives explicit extents in each link's local frame, e.g.
  link 1: `[-0.08,-0.08,-l1/2]` to `[0.08,0.08,l1/2]`, i.e. a box centered on the COM frame
  aligned with the link's local axes). For a first pass, simple cylinders along each link's
  known joint-to-joint vector (reconstructible from consecutive frame origins) are a fine,
  much simpler stand-in that still looks correct.
- **Flexible beam** (frames 5-10, i.e. 7 "node" poses once you include frame 4 as the proximal
  end): draw as a swept rectangular (or for visual simplicity, circular) cross-section between
  consecutive node frame poses (`elasticBeam.m`, Sec. above), cross-section lies in the node
  frame's local XY plane, beam axis is local Z. Default visual width/height in `elasticBeam.m`
  is `0.05` but `MBFlexibleLinkVisualization.m` (not read) likely scales this from `dBeam`;
  using the actual beam radius `2.75e-3 m` (from Sec. 2) directly as the tube radius for
  rendering is the physically correct choice and simpler than chasing the exact MATLAB default.
- Two end-effector balls (steel spheres) at beam node 1 (`r=0.03 m`) and the last beam node
  (`r=0.007 m`), per Sec. 2, render as literal spheres at those node origins, since the
  physics already treats them as point masses rigidly attached to those frames.
- World frame: z is "up" (example sets `zlim([0,1.75])` in the GitHub version's visualization
  calls; the robot stands with its base link extending from z=0 upward along the chain).

## What still needs the MATLAB reference run (not resolvable from source alone)

Everything above is traced to source and should be implementable directly. The ONE thing that
genuinely cannot be nailed down by reading alone is item 1 in Sec. 6 (gravity sign) plus general
numerical cross-validation: generate `tout`, `q`, `q_dot` (and ideally `g` for a couple of
frames) from an actual MATLAB run of this exact example and diff the Rust port's trajectory
against it at matching time steps, rather than trusting either derivation blind.
