# ELARA Optimal-Control Examples - Implementable Transcription Spec

Covers `optimal_control_planar_robot.m`, `optimal_control_rigid_robot.m`,
`optimal_control_continuum_manipulator.m`. Builds directly on `math_spec.md`'s SE3/Cayley/DEL
conventions (same `cayRTDInvSE3`, same frame/Jacobian machinery) - this document only adds what
is new for optimal control. Traced to the actual installed source; CasADi itself is NOT
installed in this MATLAB (`casadi.SX` fails to resolve), so none of this could be run live to
cross-check against a solved reference trajectory the way the simulation port was. Everything
below is read directly from source, not run.

## 1. OCPDefinition fields (full semantics)

- `tF`, `h`: horizon length and fixed step; `nSteps = round(tF/h)`, `tout = 0:h:h*nSteps`.
- `q0`, `qDot0`: initial state, enforced as **tight variable bounds** (`lb=ub=q0`) on the first
  q-column of the decision vector, not as a residual/constraint.
- `qDotF`: if non-empty, enforced through a special **final-step DEL residual**
  (`computeDELResiduumFinalStep_casadi_extKin`, mirrors the simulation's "first step" formula
  but at the other end of the horizon, folding in the desired final velocity via a discrete
  Legendre transform). If empty, the dynamics simply end at step `nSteps` and a trivial
  `u(end) - u(end-1) = 0` constraint is added instead (not exercised by any of the 3 examples;
  all three set `qDotF = zeros(nDoF,1)`).
- `qF`: if non-empty, an **explicit equality constraint** `q(:,end) == qF` is appended (comment
  in source: works better than enforcing via variable bounds). Only `optimal_control_planar_robot.m`
  sets this (`qF = [-pi/2, 0]`); the other two leave it empty.
- `uMin`/`uMax`: box bounds. Applied as **direct variable bounds** on the decision vector when
  `useSplineInputs == false` (planar robot). When `useSplineInputs == true` (rigid robot,
  continuum manipulator), the spline control points `z` are left unbounded as variables, and the
  bound is instead applied as an **explicit elementwise inequality constraint** on the
  spline-evaluated `u` at every step 1..nSteps (see `OCPConstraintFunDEL.m` lines ~119-136, and
  the final step's own copy at the end of that file). Either way the *effective* bound on `u(t)`
  is the same; which mechanism enforces it is a solver-performance detail, not a modeling one.
- `qMin`/`qMax`: always applied as direct variable bounds on every q-column.
- `wRC` (5,1) / `iRC` (5,1): running-cost weights/flags, order **[u-norm, u_dot-norm,
  u_ddot-norm, q_ddot-norm, TCP-tracking]**. `iRC` just gates which terms are summed; the weight
  itself is unconditionally computed either way (gating only matters if a term would be
  ill-defined, e.g. no TCP trajectory given).
- `wFC` (3,1) / `iFC` (3,1): final-time cost weights/flags, order **[u(end)-norm, q(end)-norm,
  TCP-final-error]**.
- `FDOrder` (2 or 4): order of the finite-difference stencil (`diffHigherOrder.m`, already in
  `elara/internal/`) used for `q_ddot` always, and for `u_dot`/`u_ddot` **only when
  `useSplineInputs == false`**; with splines, `u_dot`/`u_ddot` come from the exact spline
  derivative basis matrices instead (no finite-difference error there).
- `useSplineInputs`, `nInputSplinePoints`, `inputSplineOrder`: B-spline control parameterization,
  see Sec. 3.
- `x_TCP_F`, `R_TCP_F`: desired final TCP pose. `R_TCP_F` is declared but **never read** in
  `defineOCPSolver.m`/`OCPConstraintFunDEL.m` for any of the three examples (only `x_TCP_F`,
  the position, is used; the final-cost/constraint is a pure position error, not a pose error)
  - all three examples also leave `R_TCP_F = []` anyway, so this is moot for them, but note it's
    not wired up even if set.
- `addTCPFinalTimeConstraint`: if true, adds a **hard equality** `g_TCP_F.x - x_TCP_F == 0` in
  addition to (not instead of) whatever soft TCP cost terms are active. All three examples leave
  this `false` - TCP targeting, where used, is purely a soft cost term (JR5 and/or JF3, see Sec. 4).
- `tPreAct`, `tPostAct`, `x_TCP_waypoints`, `x_TCP_timepoints`, `x_TCP_traj`: only feed
  `generateDesiredTCPTrajLinear.m` / `OCPComputeInitialGuess_InvDyn.m` (initial-guess and
  reference-trajectory generation for the running TCP-tracking cost, used by the rigid-robot
  example). Not read in depth (out of scope per my directive) - for a port, any reasonable
  initial guess (e.g. straight-line interpolation in q from q0, or just q0 held constant) is a
  legitimate substitute; MATLAB's own initial guess only affects *convergence speed/robustness*
  of the NLP solve, never the problem definition or (for a well-posed problem) the optimum found.
- `workSpaceDef`: a `workspaceDefinition` object (box obstacles/boundaries). Default-constructed
  with `nObjects = 0`. **None of the three examples ever call `addBoxFromAxisLimits` /
  `addBoxSideLengths` on it** - confirmed by reading every line of all three scripts - so for all
  three, `dIntFun`/`dExtFun` (built from 0 boxes) produce empty constraint arrays and the
  workspace inequality constraints in `OCPConstraintFunDEL.m` are **entirely vacuous**. A port
  can skip workspace constraints completely for these three examples.
- `discretization`: `OCPIntegratorVI` (what this spec covers) or `OCPIntegratorRK` (an
  alternate RK4/RK2 direct-transcription path the examples also run once, as `OCP_ODE`, purely
  as a cross-check plotted against the VI solution - not needed for the port, and not read in
  depth here per directive scope).
- `discretization.aTrapez`: **defaults to 0.5** (true trapezoidal rule) for `OCPIntegratorVI`,
  unlike the PD simulation example which explicitly overrode it to `0`. None of the three OCP
  examples override it, so all three use `aTrapez = 0.5` - i.e. the DEL damping term is the
  symmetric/centered form here, not the one-sided form from the simulation port. Must not reuse
  the simulation's `aTrapez=0` assumption.

## 2. DEL dynamics equality constraint (per step, with control input)

Structurally identical to the simulation port's DEL residual (`math_spec.md` Sec. 4): same
per-frame inertia terms via `cayRTDInvSE3`, same gravity term, same recursive discrete-velocity
computation. Two differences for the OCP context:

1. **Control input term.** The simulation's `f_gen_k` had `-h*B(q_k)*u` with `u` always empty
   (zero contribution). Here `u_k` is a genuine free decision variable (or spline-evaluated
   value) entering through the **same** `B(q_k) = computeInputMatrixFast(g_rel_k)` mapping:
   - **Rigid/screw joints (jointType==1, the only case needed for planar-robot and
     rigid-robot):** `B(qIndex, uIndex) = 1` - i.e. `u` is *directly* the generalized torque on
     that DOF, a trivial 1:1 selection matrix. This is the only input mapping the planar and
     rigid robot examples need; fully resolved, nothing further to derive.
   - **Flexible/cable joints (jointType==2, needed only for the continuum manipulator):** a
     materially different, more involved mapping (`MBSystemNum.m` lines ~407-426): for each
     cable `iC` at frame `iFrm`, reads two "cable routing" frames `g_cm(:,:,1,iFrm,iC)` and
     `g_cm(:,:,2,iFrm,iC)` (adjacent-node cable anchor poses, presumably from
     `MBLinkCableActuationConfig.m` / `TendonActuation.m`, **not read in this pass** - flagged
     below as the one open item), computes a per-cable tangent twist via
     `cayInvSE3(g_cm_i1 \ g_rel(:,:,iFrm) * g_cm_i2) / l`, then maps that into generalized force
     via `Ba(iFrm)' * [skew(cable-anchor-midpoint) * tangent_lin; tangent_lin]`, scaled by
     `-l/norm(tangent_lin)`. This is a real lever-arm/moment-arm computation for a cable
     routed near (not necessarily through) the beam's centerline. **Do not attempt the
     continuum-manipulator port without first reading `MBLinkCableActuationConfig.m`,
     `TendonActuation.m`, and the cable-routing setup inside `systemDef_continuum_manipulator.m`
     itself (not yet read) to get the actual `g_cm` values this example uses.**

2. **Scale/sign convention.** The OCP's own residual functions
   (`computeDELResiduum_extKin.m`, `computeDELResiduumFirstStep_extKin.m`) use a *different*
   internal scaling than `integrateMBSDynamics_Broyden.m` (e.g. the general-step version has NO
   `h` multiplying the stiffness/input term and an overall `-1` sign flip at the end; the
   first-step version matches the simulation's `h*f_gen - p_0` form with no final sign flip).
   **This does not matter for a port**: `constraint == 0` is invariant to an overall nonzero
   scale or sign, so reusing the already-validated simulation-port residual formula (from
   `src/dynamics.rs`), just adding the `-h*B(q_k)*u_k` input term that the simulation path always
   zeroed out, is a safe and simpler choice than replicating MATLAB's particular scaling choice
   per formula. Keep the convention *internally consistent* (same sign/scale for first/general/
   final step within your own port), which is simpler to guarantee by reusing one formula
   throughout than by copying three differently-scaled MATLAB formulas.

Per-step dimension: the DEL residual is `nDoF`-dimensional equality constraint at every interior
step (`k = 1 .. nSteps-1`, i.e. `nSteps-1` vector constraints), plus one more at the first step
and one more at the final step (using the `qDotF`-seeded final-step formula) - total
`nSteps+1` vector equality constraints of size `nDoF`, matching one DEL constraint per
state-to-state transition in the horizon.

## 3. Decision variables

**Without spline parameterization** (`optimal_control_planar_robot.m`): decision vector is
simply `[q_0; u_0; q_1; u_1; ...; q_N; u_N]` (or any other fixed bijection with the same
content) - `q` and `u` are independent per-step decision variables, `nDoF + nInputs` variables
per step node, `(nSteps+1)*(nDoF+nInputs)` total. `MATLAB`'s actual `getXVecIndicesSpline`-driven
interleaving of q/u (or q/spline-point) slots inside the flat NLP vector is purely a
sparsity/performance choice for CasADi/IPOPT - **the variable ordering has zero effect on the
solved trajectory**, only on solver internals. A port is free to use the simplest possible
layout (e.g. all q's, then all u's, or one struct-of-arrays).

**With spline parameterization** (rigid robot, continuum manipulator): `u` is not a free
per-step variable; instead `z` (shape `nInputs x nInputSplinePoints`) are the free decision
variables (control points), and `u(t) = B(t) * z` via the clamped B-spline basis matrix (Sec. 3a
below). Decision vector = `[q_0..q_N, z_1..z_{nInputSplinePoints}]` in any order -
`nDoF*(nSteps+1) + nInputs*nInputSplinePoints` variables total.

### 3a. B-spline basis (`computeBSplineBasisMatrix.m`)

Standard **clamped** (open) B-spline, Cox-de Boor recursion, degree `p = inputSplineOrder`
(default/always 3, cubic, in all three examples), `nCP = nInputSplinePoints` control points.
Knot vector: `(p+1)`-fold knot at both `t0` and `tEnd`, uniformly spaced interior knots:
```
K = nCP - 1;  m = K + p + 1;  nT = m + 1
T = [repeat(t0, p), linspace(t0, tEnd, nT - 2p), repeat(tEnd, p)]
```
Zero-degree basis `B0(k,i) = 1` if `t_k` is in `[T(i), T(i+1))`, with a special case to make the
*last* basis function active at exactly `t = tEnd` (multiple end knots would otherwise make it
half-open on the wrong side there). Then the standard Cox-de Boor recursion up to degree `p`,
carrying along the closed-form first/second derivative recursions at the same time (this is a
completely standard algorithm - see *The NURBS Book* ch. 2.1 - safe to implement directly in
Rust with no external crate; the whole function is ~70 lines of pure arithmetic, no linear
solves). `u(t) = B(t) @ z`, `u_dot(t) = B_dt(t) @ z`, `u_ddot(t) = B_ddt(t) @ z` - these are used
**directly** (not through any finite-difference) everywhere a spline is active.

## 4. Objective function (exact, from `defineOCPSolver.m` lines ~164-219)

Let `trapz(x)` denote `OCP.h * (sumsqr(x[:,0])/2 + sumsqr(x[:,1:-1]) + sumsqr(x[:,-1])/2)` - a
plain trapezoidal-rule quadrature of the columnwise sum-of-squares (`sumsqr` sums the squares of
*all* elements of its argument, i.e. this is `h * integral(||x(t)||^2 dt)` via the trapezoidal
rule, no `1/2` factor anywhere in `trapz` itself).

**Running cost** (`JR`, 5 terms, `wR`=`OCP.wRC`, each gated by `OCP.iRC(i)`):
```
JR[1] = wR[1] * trapz(u)          # u as used (spline-evaluated or raw decision var)
JR[2] = wR[2] * trapz(u_dot)      # spline derivative basis, or FD of order FDOrder if no spline
JR[3] = wR[3] * trapz(u_ddot)     # likewise
JR[4] = wR[4] * trapz(q_ddot)     # ALWAYS via finite difference (q is never spline-parameterized)
JR[5] = wR[5] * sum_{k=0}^{nSteps} fTrap(k) * h * 0.5 * ||x_TCP_traj[:,k] - TCP_position(q_k)||^2
        where fTrap(k) = 0.5 if k==0 or k==nSteps, else 1   (hand-rolled trapezoidal weight)
```
**Note the asymmetry**: JR[1..4] have no extra `1/2` (it's baked into nothing - `trapz` itself
has none either); JR[5] has an explicit extra `* 0.5` that `trapz` does NOT supply (JR[5] is
computed by a bespoke loop, not via `OCP.discretization.integrateCostFunctionValue`). Easy to
get wrong by analogy - implement JR[5] with its own `0.5` factor, JR[1..4] without.

**Final-time cost** (`JF`, 3 terms, `wF`=`OCP.wFC`, each gated by `OCP.iFC(i)`), all WITH an
explicit `0.5`:
```
JF[1] = wF[1] * 0.5 * ||u(:,end)||^2
JF[2] = wF[2] * 0.5 * ||q(:,end)||^2        # raw norm of q(T), not q(T)-qF - a regularizer, not
                                              # a tracking error (only matters when OCP.qF is NOT
                                              # also hard-constrained; none of the 3 examples use
                                              # iFC[2], so this is moot for them, but implement
                                              # literally if porting the field)
JF[3] = wF[3] * 0.5 * ||x_TCP_F - TCP_position(q(:,end))||^2
```
`J = sum(JR[i] for i where iRC[i]) + sum(JF[i] for i where iFC[i])`.

`TCP_position(q)` = translation part of `g(indexTCPFrame) * g_B_TCP` (forward kinematics of the
TCP frame, composed with the fixed TCP offset already in the system definition) - exactly the
same quantity the simulation port's rendering layer already computes per frame, just picked out
for one specific frame index.

Per-example `iRC`/`iFC` (read directly from the three scripts, for a quick sanity check when
porting):
- **planar_robot**: `wRC = [1/2, 0, 0, 0, 0]`, `iRC = logical(wRC)` -> only `JR[1]` (pure
  control-effort minimization). `wFC = iFC = zeros(3,1)` -> no final cost at all. Dynamics
  constraint is DEL + the hard `qF` equality + the `qDotF`-seeded final-step DEL. This is the
  cleanest, smallest problem of the three (good first port target): 2 DOF, no TCP cost, no
  splines, point-to-point with a hard final-state constraint.
- **rigid_robot**: `wRC = [1e-2/2, 0, 0, 0, 1e3]` -> control effort + TCP-trajectory-tracking
  running cost (needs `x_TCP_traj`, generated by `generateDesiredTCPTrajLinear.m`, not read in
  depth - for a port, any reasonable straight-line-in-task-space interpolation from the initial
  to final TCP position over `[tPreAct, tF-tPostAct]`, held constant outside that window, is a
  faithful reimplementation of what that function almost certainly does given its name and the
  `tPreAct`/`tPostAct` fields it consumes; worth a direct read before implementing if exactness
  matters). `iFC = zeros(3,1)` -> no final cost. No `qF` set -> no hard final-configuration
  constraint, only the `qDotF`-seeded final-step DEL.
- **continuum_manipulator**: `wRC = [1e-3, 1e-3, 1e-3, 1e-2, 0]` -> effort + smoothness
  (u_dot, u_ddot) + q_ddot penalty, no running TCP cost. `wFC = [0, 0, 1e-3*1e8]`, `iFC` only
  term 3 active -> pure **soft final-TCP-position cost** (no hard constraint, no `qF`). This is
  the one that also needs the cable-actuation input mapping from Sec. 2.

## 5. Constraint summary per example

| | planar_robot | rigid_robot | continuum_manipulator |
|---|---|---|---|
| DEL dynamics (eq., every step) | yes | yes | yes (needs cable `B` matrix) |
| Final-step DEL w/ `qDotF` | yes | yes | yes |
| Hard final config `q(T)=qF` | **yes** | no | no |
| Hard final TCP (`addTCPFinalTimeConstraint`) | no | no | no |
| Soft TCP running cost (JR5) | no | yes | no |
| Soft TCP final cost (JF3) | no | no | yes |
| `u` bounds | var. bounds (no spline) | nonlinear constraint (spline) | nonlinear constraint (spline) |
| `q` bounds | yes (`[-inf,-4]`..`[inf,0.1]`-ish, see script) | yes (`±2pi`) | not set in script (unbounded) |
| Workspace obstacles | inactive (0 boxes) | inactive (0 boxes) | inactive (0 boxes) |

## 6. New robot system definitions (exact numbers)

Both probed live via `matlab -batch` against the installed toolbox (same method as the
simulation port) - see `reference/ocp_planar_frame_data.json` and
`reference/ocp_rigid_frame_data.json` for the full per-frame `MGen`/`X`/`g_ref` dumps, ready to
feed the same JSON->Rust-constants conversion script already used for the simulation port
(`reference/gen_rust_data.py` pattern).

**Planar 2-link pendulum** (`systemDef_planarNLinkPendulum("nLinks", 2, "d", 0)`): 2 identical
rigid links, **unit mass each** (`m=1`, by construction: radius solved from
`r = sqrt(m/(rho*pi*l))` with `rho=2700`, `l=1`), joint axis **local Y** (planar motion in the
X-Z plane, not Z like the other two robots), serial chain, `d=0` (frictionless, as passed by the
OCP example). Confirmed via probe: `m = [1, 1]`, `jointType = [1,1]`, screw axes
`X = [[0,1,0,0,0,-0.5], [0,1,0,0,0,-0.5]]` (both links: pure Y-rotation, translation component
`-0.5` along the twist's linear part = `-l/2`, consistent with `g_J_B = g_J1_COM =
SE3(I,[l/2,0,0])` and the AdInv transform collapsing that offset onto the Y-axis lever arm).
`nDoF=2`, `nInputs=2` (both joints actuated), `indexTCPFrame=0` -> **no TCP frame explicitly
defined**, so `defineOCPSolver.m` falls back to the *last* frame (frame 2) with a warning; not
that it matters, since neither the planar example's cost nor constraints use a TCP term.
Inertia: `I_sym = 1/2*r^2` (about the link's own long axis, here the out-of-plane rotation axis
doesn't apply the same way as the other robots since the joint axis IS the symmetry-breaking
axis - just use the two diag values directly from the JSON rather than re-deriving which is
which). `g0 = I` (world frame = base frame).

**Rigid 3-link lab robot** (`systemDef_rigid_robot("d", 0)`): structurally identical in *form*
to links 1-3 of the rigid-flexible robot already ported (same density `rho=2700`, same radii
`0.05/0.03/0.03`, same inertia formulas, same `g_J_const` 90-degree twist between link 1 and 2,
same joint-axis-Z convention) but with **different link lengths**: `l = [0.2, 0.5, 0.5]` (vs.
`[0.25, 0.4, 0.4]` in the rigid-flexible robot) and `d=0` (vs. `0.01` there, since the OCP
example explicitly zeroes joint friction for a cleaner optimal-control problem). Probed masses:
`m = [4.2412, 3.8170, 3.8170]` kg. `nDoF=3`, `nInputs=3`, `indexTCPFrame=3` (last link has
`hasTCP=true`), `g_B_TCP = SE3(I, [0.25,0,0])` = `l(3)/2` forward of link 3's own frame origin
(i.e. the TCP sits at the *tip* of link 3, which is `+l(3)/2` ahead of the COM-centered local
frame). Screw axes (probed): `X = [[0,0,1,0,0,0], [0,0,1,0,0.25,0], [0,0,1,0,0.25,0]]` - link 1 is
pure Z-rotation at the origin (no lever-arm component, matching how link 1 of the rigid-flexible
robot worked too); links 2 and 3 both carry a `0.25` linear-twist component (`= l(2)/2` /
`l(3)/2`? note both show `0.25` even though `l(2)=l(3)=0.5`, consistent with `l/2`) in the same
pattern as the rigid-flexible robot's links 2/3. This system has **no flexible component at all**
- it is exactly a 3-DOF open rigid chain, nothing new needed beyond what `src/kinematics.rs`
already implements for rigid/screw frames (frames 1-3 of the existing rigid-flexible port are
literally this same formula family).

**Continuum manipulator** (`systemDef_continuum_manipulator("nSeg", 4)`): **not yet read** -
out of scope for this pass per the time I had; needed before that example can be ported at all,
both for its own link/beam parameters AND for the cable-actuation `g_cm` routing data Sec. 2
flags as the real open item. Recommend reading `systemDef_continuum_manipulator.m`,
`MBLinkCableActuationConfig.m`, and `TendonActuation.m` together in one pass before attempting
this third example - it is a strictly larger lift than planar_robot or rigid_robot.

## 7. What I could not pin down / did not attempt (be explicit about these before porting)

1. **Continuum manipulator's cable actuation** (Sec. 2, item 2) - the `g_cm` cable-routing
   frames are not yet sourced. Do planar_robot and rigid_robot first; they need none of this.
2. **`generateDesiredTCPTrajLinear.m`** (rigid-robot's reference TCP trajectory for JR5) - not
   read; inferred likely behavior from field names only (see Sec. 4). Read it directly before
   trusting the inference if exact cost-value matching to MATLAB is wanted.
3. **`OCPComputeInitialGuess_InvDyn.m` / `computeOptimalSteadyStateInputsTCPPos.m`** - not read
   in depth (affects only initial-guess quality / solver convergence speed, never the problem
   definition or, for a well-posed problem, the optimum itself). Any reasonable initial guess is
   fine for a from-scratch Rust solver.
4. **`OCPIntegratorRK.m` / `RKStepOCP.m`** - the alternate ODE/RK4 transcription used for the
   examples' secondary cross-check solve (`OCP_ODE`). Not needed for the port; not read.
5. **No live numerical cross-check was possible** for this spec (CasADi is not installed in this
   MATLAB, so none of the three `optimal_control_*.m` scripts can actually be run end-to-end
   here) - unlike the simulation port, there is no MATLAB-solved reference trajectory to validate
   a Rust port against yet. Getting CasADi+IPOPT installed and actually running these scripts to
   capture a real solved `(q_sol, u_sol)` trajectory (the same way `vi_q_qdot.csv` was captured
   for the simulation) should happen before trusting a Rust NLP solver's output as "matching
   MATLAB" - right now there is nothing concrete to match against, only the transcription itself
   traced faithfully to source.
6. **No NLP solver is specified here** - this document covers the math *problem* (cost,
   constraints, variables), not how to solve it. Porting IPOPT itself is not realistic in Rust/
   WASM (it depends on MUMPS/HSL sparse linear algebra); a custom equality-constrained
   Newton/SQP or interior-point method using exact (e.g. forward-mode dual-number) derivatives of
   the DEL residual and cost is the realistic path, and is a separate, substantial design
   decision left to whoever implements this next.
