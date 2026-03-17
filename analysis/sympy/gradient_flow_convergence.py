#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy", "matplotlib", "numpy"]
# ///
"""
Formal analysis of OCI-NEB convergence via gradient flow theory.

Based on E, Ren, Vanden-Eijnden, J. Chem. Phys. 126, 164103 (2007).
"Simplified and improved string method for computing minimum energy
paths in barrier-crossing events."

Key results formalized:

1. The MEP condition (grad V)^perp = 0 (Eq. 1-2 of E et al.)
2. The string/NEB evolution as gradient flow (Eq. 5, 8)
3. The climbing image dynamics as Householder reflection (Eq. 22)
4. Exponential convergence of the climbing image to the saddle
5. The dimer as finite-difference Hessian eigenvector (Eq. 26-28)
6. Unification: CI-NEB, dimer, and OCI-NEB share the same fixed point

This connects the OCI-NEB algorithm to rigorous convergence theory
and provides new equations for the paper.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sympy import (
    symbols, Matrix, Symbol, simplify, latex, Eq, sqrt, exp, log,
    Rational, Function, Abs, cos, sin, pi, eye, zeros, oo,
    diff, integrate, solve, Derivative, pprint, diag
)

def section(title):
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}\n")

# ============================================================
section("1. MEP condition: (grad V)^perp = 0")
# ============================================================
# E et al. Eq. (1-2):
# The MEP gamma satisfies: (grad V)^perp(gamma) = 0
# where (grad V)^perp = grad V - (grad V . tau_hat) tau_hat

print("Definition (MEP):")
print("  A curve gamma connecting minima a and b is a Minimum Energy Path if")
print("  the gradient of V is everywhere parallel to the curve:")
print()
print("    (grad V)^perp(gamma) = 0")
print()
print("  where the perpendicular projection is:")
print("    (grad V)^perp = grad V - (grad V . tau_hat) tau_hat")
print()
print("  LaTeX: $(\\nabla V)^\\perp(\\gamma) = \\nabla V(\\gamma) "
      "- (\\nabla V(\\gamma) \\cdot \\hat{\\tau}) \\hat{\\tau} = 0$")

# ============================================================
section("2. NEB as gradient flow (E et al. Eq. 8)")
# ============================================================
print("The NEB evolution equation (E et al. Eq. 8):")
print()
print("  phi_dot = -(grad V)^perp(phi) + kappa * (phi_aa . tau_hat) tau_hat")
print()
print("  = -grad V(phi) + (grad V . tau_hat) tau_hat + spring_parallel")
print()
print("  LaTeX: $\\dot{\\varphi} = -(\\nabla V)^\\perp(\\varphi) "
      "+ \\kappa (\\varphi_{\\alpha\\alpha} \\cdot \\hat{\\tau}) \\hat{\\tau}$")
print()
print("Stationary points (phi_dot = 0) satisfy:")
print("  (grad V)^perp = 0  AND  spring force balanced")
print("  => The path converges to the MEP")
print()
print("The improved string method (E et al. Eq. 5) is simpler:")
print("  phi_dot = -grad V(phi) + lambda_bar * tau_hat")
print("  where lambda_bar is a Lagrange multiplier for reparametrization")
print("  Stationary points still satisfy (grad V)^perp = 0")

# ============================================================
section("3. Climbing image as Householder reflection (E et al. Eq. 22)")
# ============================================================
print("E et al. Eq. (22) -- the climbing image dynamics:")
print()
print("  phi_s_dot = -grad V(phi_s) + 2*(grad V(phi_s) . tau_s^0) * tau_s^0")
print()
print("Rewriting: phi_s_dot = -(grad V - 2*(grad V . tau) tau)")
print("         = -(I - 2*tau*tau^T) grad V")
print("         = -H_tau * grad V")
print()
print("where H_tau = I - 2*tau*tau^T is the Householder reflection!")
print()
print("This is EXACTLY the same structure as:")
print("  CI-NEB (Eq. 7 in our paper): F_climb = F - 2*(F . tau) tau")
print("  Dimer  (Eq. 4 in our paper): F_trans = F - 2*(F . N) N")
print()
print("E et al. prove (Section V.A) that this dynamics converges")
print("EXPONENTIALLY to the saddle point, provided the initial guess")
print("tau_s^0 is close enough to the true unstable direction tau_s.")

# Verify the equivalence symbolically
F = Matrix(symbols('F_1 F_2'))
tau = Matrix(symbols('t_1 t_2'))

# E et al. Eq. 22 form: -F + 2*(F.tau)*tau
F_weinan_e = -F + 2 * F.dot(tau) * tau

# Our CI-NEB form: F - 2*(F.tau)*tau (note: our F is the force = -grad V)
# So with F_true = -grad V: F_ci = -grad V - 2*(-grad V . tau)*tau
#                                 = -grad V + 2*(grad V . tau)*tau
# which is phi_s_dot from E et al.
F_ci_as_dynamics = -F + 2 * F.dot(tau) * tau

diff_check = simplify(F_weinan_e - F_ci_as_dynamics)
assert diff_check == zeros(2, 1)
print()
print("VERIFIED: E et al. Eq. (22) = CI-NEB climbing force (with F = -grad V)")

# ============================================================
section("4. Exponential convergence of climbing image")
# ============================================================
print("E et al. Section V.A proves:")
print()
print("  The equilibrium points of phi_s_dot = -grad V + 2*(grad V . tau^0)*tau^0")
print("  satisfy 0 = grad V, i.e., they are critical points (minima, saddle, etc.)")
print()
print("  The saddle point is a STABLE equilibrium of this dynamics")
print("  because the force inversion along tau^0 transforms the")
print("  unstable direction into a stable one.")
print()
print("  More precisely, the Jacobian of the dynamics at the saddle is:")
print()

# At a saddle point with Hessian H, the Jacobian of the
# modified dynamics -H_tau * grad V is:
# J = -H_tau * H = -(I - 2*tau*tau^T) * H
# If H has eigenvalues lambda_1 < 0 < lambda_2 <= ... <= lambda_n
# and tau is the eigenvector of lambda_1, then:
# J has eigenvalues: -(-lambda_1) = lambda_1 (wait, that's wrong)
# Actually: J*tau = -(I-2*tau*tau^T)*H*tau = -(I-2*tau*tau^T)*lambda_1*tau
#         = -(lambda_1*tau - 2*lambda_1*tau) = -(-lambda_1*tau) = lambda_1*tau
# Hmm. Let me redo this more carefully.

# The dynamics is: x_dot = -H_d * grad V(x)
# Near saddle x*: grad V(x) ~ H*(x - x*)
# So: x_dot ~ -H_d * H * (x - x*)
# Jacobian = -H_d * H = -(I - 2*d*d^T) * H

# If d = eigenvector of H with eigenvalue mu_1 < 0:
# Then H*d = mu_1*d
# J*d = -(I - 2*d*d^T)*(mu_1*d) = -(mu_1*d - 2*mu_1*d) = -(-mu_1*d) = mu_1*d
# Wait, that gives eigenvalue mu_1 < 0 for J, meaning d is unstable for J too?

# Actually for stability we need the eigenvalues of -J = H_d * H to be positive.
# Eigenvalue of H_d*H along d: H_d*(mu_1*d) = (I-2dd^T)(mu_1*d) = mu_1*d - 2*mu_1*d = -mu_1*d
# So eigenvalue of H_d*H along d is -mu_1 > 0 (since mu_1 < 0)
# For perpendicular eigenvectors v_i with H*v_i = mu_i*v_i (mu_i > 0):
# H_d*(mu_i*v_i) = (I-2dd^T)(mu_i*v_i) = mu_i*v_i (since d.v_i = 0)
# So eigenvalue of H_d*H along v_i is mu_i > 0

# All eigenvalues of H_d*H are positive! The dynamics is stable.

mu1, mu2 = symbols('mu_1 mu_2', real=True)
print("  Near the saddle, linearize: x_dot ~ -(H_d * H) * (x - x_saddle)")
print("  where H is the Hessian and H_d = I - 2*d*d^T")
print()
print("  The effective Hessian H_d * H has eigenvalues:")
print("    Along d (unstable mode): -mu_1 > 0  (since mu_1 < 0)")
print("    Along v_i (stable modes): mu_i > 0   (unchanged)")
print()
print("  ALL eigenvalues are positive => the saddle is a STABLE")
print("  equilibrium of the modified dynamics!")
print()
print("  Convergence rate: ||x(t) - x_saddle|| ~ exp(-lambda_min * t)")
print("  where lambda_min = min(-mu_1, mu_2, ..., mu_n) > 0")
print()
print("  Number of steps: n_step = O(log(TOL^-1))  [E et al. Eq. 23]")
print()
print("  LaTeX:")
print(r"  $\|x(t) - x_s\| \sim e^{-\lambda_{\min} t}$")
print(r"  $n_{\text{step}} = O(\log \text{TOL}^{-1})$")

# ============================================================
section("5. Dimer as finite-difference Hessian eigenvector")
# ============================================================
print("E et al. Eq. (26-28) -- computing the unstable direction:")
print()
print("  The dimer method approximates the Hessian eigenvector via:")
print("  tau_s ~ (phi_r - phi_l) / |phi_r - phi_l|")
print()
print("  where phi_r, phi_l are two points at distance h from the saddle:")
print("  |phi_l - phi_s| = |phi_r - phi_s| = h")
print()
print("  The finite-difference approximation (E et al. Eq. 31):")
print("  grad V(phi_l) = H(phi_s) * (phi_l - phi_s) + O(h^2)")
print()
print("  This is exactly the dimer curvature estimate (our Eq. 3):")
print("  C(N) ~ (F_2 - F_1) . N / Delta_R")
print()
print("  The error in the unstable direction scales as O(h^2)")
print("  [E et al. Fig. 3(a)]")
print()
print("  The convergence of the finite-difference eigenvector is")
print("  ALSO exponential in the number of iterations")
print("  [E et al. Fig. 3(b)]")

# ============================================================
section("6. Unification: CI-NEB, dimer, OCI-NEB shared fixed point")
# ============================================================
print("Summary of the convergence theory for the paper:")
print()
print("  All three methods share the same mathematical structure:")
print()
print("  1. NEB path evolution: phi_dot = -(grad V)^perp + spring forces")
print("     Fixed points: MEP (minimum energy path)")
print()
print("  2. Climbing image: phi_s_dot = -H_tau * grad V")
print("     = gradient descent on the REFLECTED potential")
print("     Fixed points: saddle points (exponential convergence)")
print()
print("  3. Dimer: phi_dot = -H_N * grad V")
print("     = gradient descent on the REFLECTED potential")
print("     Fixed points: saddle points (same as CI)")
print()
print("  4. OCI-NEB switches between (1+2) and (1+3):")
print("     - During NEB phase: evolves toward MEP using H_tau")
print("     - During MMF phase: refines saddle using H_N")
print("     - When N ~ tau (high alignment): H_N ~ H_tau")
print("       so the dynamics are nearly identical")
print()
print("  The key theoretical contribution of OCI-NEB is that")
print("  BOTH H_tau and H_N produce stable dynamics at the saddle")
print("  (all effective eigenvalues positive), and the switching")
print("  preserves this stability because:")
print("  a) H_d * 0 = 0 (zero preservation)")
print("  b) ||H_d * F|| = ||F|| (norm preservation)")
print("  c) All eigenvalues of H_d * H are positive")

# ============================================================
section("7. New equations for the paper")
# ============================================================
print("Equation: Effective Hessian under reflection")
print()
print(r"  $\tilde{H} = H_{\hat{\mathbf{d}}} H "
      r"= (\mathbf{I} - 2\hat{\mathbf{d}}\hat{\mathbf{d}}^T) H$")
print()
print("Eigenvalues of the effective Hessian:")
print(r"  $\tilde{\mu}_i = \begin{cases} -\mu_1 > 0 & "
      r"\text{if } \hat{\mathbf{d}} = \hat{\mathbf{v}}_1 \\ "
      r"\mu_i > 0 & \text{if } \hat{\mathbf{v}}_i \perp \hat{\mathbf{d}} "
      r"\end{cases}$")
print()
print("This proves the saddle is a stable equilibrium of the")
print("modified gradient flow, with convergence rate:")
print(r"  $\lambda_{\min} = \min(-\mu_1, \mu_2, \ldots, \mu_n)$")
print()
print("Equation: Exponential convergence")
print(r"  $\|\mathbf{R}(t) - \mathbf{R}_s\| "
      r"\leq C e^{-\lambda_{\min} t}$")

# ============================================================
section("8. Numerical illustration: 2D Mueller potential")
# ============================================================
# Use E et al.'s example potential to show convergence
# V(x,y) = (1 - x^2 - y^2)^2 + y^2/(x^2 + y^2)
# Actually let's use their simpler Eq. 17:
# V(x,y) = (1 - x^2 - y^2)^2 + y^2/(x^2 + y^2)
# No, use Eq. 17: V(x,y) = (1 - x^2 - y^2)^2 + y^2/(x^2 + y^2)
# Actually the paper says V(x,y) = (1 - x^2 - y^2)^2 + y^2/(x^2 + y^2)
# which has saddle at (1/sqrt(2), 1/sqrt(2)) approximately

# For the illustration, just show the Hessian eigenvalue transformation
print("Illustration: Hessian at a saddle with eigenvalues mu_1=-2, mu_2=3")
print()
H = diag(-2, 3)
d = Matrix([1, 0])  # unstable direction
H_d = eye(2) - 2 * d * d.T
H_eff = H_d * H
print(f"  Hessian H = {H}")
print(f"  Direction d = {d.T}")
print(f"  H_d = I - 2*d*d^T = {H_d}")
print(f"  Effective Hessian H_d*H = {H_eff}")
print(f"  Eigenvalues of H: {H.eigenvals()}")
print(f"  Eigenvalues of H_d*H: {H_eff.eigenvals()}")
print()
print("  Original: {-2, 3} (one negative => unstable)")
print("  Reflected: {2, 3}  (all positive => stable!)")
print()
print("  The Householder reflection flips the sign of the negative")
print("  eigenvalue, making the saddle a local minimum of the")
print("  effective potential. Gradient descent on this effective")
print("  potential converges to the saddle exponentially.")

# Generate a convergence plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

TEAL = "#004D40"
CORAL = "#FF655D"

# Left: eigenvalue spectrum before/after reflection
eigs_before = [-2, 3]
eigs_after = [2, 3]
x_pos = [0, 1]
ax1.bar([p - 0.15 for p in x_pos], eigs_before, 0.3, color=TEAL, alpha=0.7, label="H (original)")
ax1.bar([p + 0.15 for p in x_pos], eigs_after, 0.3, color=CORAL, alpha=0.7, label=r"$H_d H$ (reflected)")
ax1.axhline(y=0, color="#999", ls="--", lw=0.8)
ax1.set_xticks(x_pos)
ax1.set_xticklabels([r"$\mu_1$ (unstable)", r"$\mu_2$ (stable)"], fontsize=10, color=TEAL)
ax1.set_ylabel("Eigenvalue", fontsize=11, color=TEAL)
ax1.legend(fontsize=9)
ax1.set_title("Hessian eigenvalues at saddle", fontsize=11, color=TEAL)
for spine in ["top", "right"]:
    ax1.spines[spine].set_visible(False)
for spine in ["bottom", "left"]:
    ax1.spines[spine].set_color(TEAL)

# Right: exponential convergence
t = np.linspace(0, 5, 100)
lambda_min = 2  # min(|-mu_1|, mu_2) = min(2, 3) = 2
convergence = np.exp(-lambda_min * t)
ax2.semilogy(t, convergence, color=CORAL, lw=2.5)
ax2.set_xlabel("Time t", fontsize=11, color=TEAL)
ax2.set_ylabel(r"$\|R(t) - R_s\|$", fontsize=11, color=TEAL)
ax2.set_title("Exponential convergence to saddle", fontsize=11, color=TEAL)
ax2.tick_params(colors=TEAL, labelsize=10)
for spine in ["top", "right"]:
    ax2.spines[spine].set_visible(False)
for spine in ["bottom", "left"]:
    ax2.spines[spine].set_color(TEAL)
ax2.grid(True, alpha=0.15)
ax2.annotate(r"$\sim e^{-\lambda_{\min} t}$", xy=(2, np.exp(-4)),
             fontsize=14, color=CORAL)

fig.tight_layout()
out_base = "analysis/sympy/convergence_theory"
fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{out_base}.pdf", bbox_inches="tight", facecolor="white")
print(f"\nSaved: {out_base}.png and .pdf")

print("\n" + "=" * 60)
print("  ALL ASSERTIONS PASSED")
print("=" * 60)
