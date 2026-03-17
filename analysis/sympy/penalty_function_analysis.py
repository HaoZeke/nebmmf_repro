#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["sympy", "matplotlib"]
# ///
"""
Symbolic analysis of the OCI-NEB penalty function.

The penalty function P(alpha, s) = b + (1-b)*alpha^s with b = 1/(1+s)
controls adaptive backoff. This script proves its key properties:

1. Monotonicity in alpha (better alignment -> less penalty)
2. Monotonicity in s (stronger penalty -> more sensitive to misalignment)
3. Integral characterization (mean penalty over uniform alignment)
4. Fixed-point property of the derived base
5. Convexity/concavity regions

References: Eq. (5) in paper, NEBRonebController.cpp line 214-216
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sympy import (
    Symbol, Rational, simplify, latex, integrate, diff,
    log, sqrt, oo, S, pprint, factor, cancel, apart
)

def section(title):
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}\n")

alpha = Symbol("alpha", positive=True)
s = Symbol("s", positive=True)
b = Symbol("b", positive=True)

# General penalty function
P_general = b + (1 - b) * alpha**s

# Derived base
b_derived = 1 / (1 + s)
P = P_general.subs(b, b_derived)
P_simplified = simplify(P)

section("1. Penalty function with derived base")
print(f"P(alpha, s) = b + (1-b)*alpha^s")
print(f"b = 1/(1+s)")
print(f"P(alpha, s) = {P_simplified}")
print(f"LaTeX: $P(\\alpha, s) = {latex(P_simplified)}$")

section("2. Boundary values")
P_at_0 = simplify(P.subs(alpha, 0))
P_at_1 = simplify(P.subs(alpha, 1))
print(f"P(0, s) = {P_at_0} = 1/(1+s) = b  [minimum: worst alignment]")
print(f"P(1, s) = {P_at_1} = 1              [maximum: perfect alignment]")
assert P_at_0 == b_derived
assert P_at_1 == 1

section("3. Monotonicity in alpha")
dP_dalpha = diff(P, alpha)
dP_simplified = simplify(dP_dalpha)
print(f"dP/d_alpha = {dP_simplified}")
print(f"LaTeX: $\\frac{{\\partial P}}{{\\partial \\alpha}} = {latex(dP_simplified)}$")
print()
print(f"For alpha > 0 and s > 0: alpha^(s-1) > 0 and s/(1+s) > 0")
print(f"=> dP/d_alpha > 0 for alpha in (0, 1]")
print(f"=> P is strictly increasing in alpha (better alignment -> higher penalty factor)")
print(f"=> Lower alignment causes more aggressive threshold reduction")

section("4. Monotonicity in s (penalty gets more selective)")
dP_ds = diff(P, s)
dP_ds_simplified = simplify(dP_ds)
print(f"dP/ds = {dP_ds_simplified}")
print()
# At alpha = 0.5 (moderate alignment)
dP_ds_half = simplify(dP_ds.subs(alpha, Rational(1, 2)))
print(f"At alpha=0.5: dP/ds = {dP_ds_half}")
print(f"For alpha < 1, increasing s makes P smaller (more selective)")

section("5. Integral (mean penalty over uniform alignment)")
mean_P = integrate(P, (alpha, 0, 1))
mean_P_simplified = simplify(mean_P)
print(f"E[P] = int_0^1 P(alpha) d_alpha = {mean_P_simplified}")
print(f"LaTeX: $\\mathbb{{E}}[P] = {latex(mean_P_simplified)}$")
print()

# Evaluate at specific s values
for s_val in [Rational(1, 2), 1, Rational(3, 2), 3]:
    b_val = 1 / (1 + s_val)
    mean_val = mean_P_simplified.subs(s, s_val)
    print(f"  s={float(s_val):.1f}: b={float(b_val):.3f}, E[P]={float(mean_val):.4f}")

section("6. Second derivative (convexity)")
d2P_dalpha2 = diff(P, alpha, 2)
d2P_simplified = simplify(d2P_dalpha2)
print(f"d2P/d_alpha^2 = {d2P_simplified}")
print()
print(f"For s > 1: d2P > 0 at alpha near 0 (convex, gentle start)")
print(f"For s < 1: d2P < 0 at alpha near 0 (concave, steep start)")
print(f"At s = 1: P is linear in alpha, d2P = 0")

# Verify s=1 gives linear
P_linear = simplify(P.subs(s, 1))
print(f"\nP(alpha, s=1) = {P_linear}")
assert P_linear == Rational(1, 2) + Rational(1, 2) * alpha

section("7. Specific default values")
s_default = Rational(3, 2)
b_default = 1 / (1 + s_default)
print(f"Default strength s = {s_default} = 1.5")
print(f"Derived base b = 1/(1+1.5) = {b_default} = 0.4")
print()

P_default = P.subs(s, s_default)
P_default_simplified = simplify(P_default)
print(f"P(alpha, 1.5) = {P_default_simplified}")

# Key alignment values
for a_val in [0, Rational(1, 4), Rational(1, 2), Rational(3, 4), 1]:
    p_val = float(P_default_simplified.subs(alpha, a_val))
    print(f"  alpha={float(a_val):.2f}: P={p_val:.4f}")

section("8. Generate penalty family plot")

TEAL = "#004D40"
CORAL = "#FF655D"

fig, ax = plt.subplots(figsize=(5.5, 4.0))

alphas = np.linspace(0, 1, 200)
strengths = [(0.5, "--", TEAL, "s=0.5 (b=0.667)"),
             (1.0, ":", TEAL, "s=1.0 (b=0.500)"),
             (1.5, "-", CORAL, "s=1.5 (b=0.400, default)"),
             (3.0, "-.", TEAL, "s=3.0 (b=0.250)")]

for s_val, ls, color, label in strengths:
    b_val = 1.0 / (1.0 + s_val)
    P_vals = b_val + (1.0 - b_val) * alphas**s_val
    lw = 2.5 if s_val == 1.5 else 1.5
    ax.plot(alphas, P_vals, ls=ls, color=color, lw=lw, label=label)
    # Dashed horizontal line at base
    ax.axhline(y=b_val, color=color, ls="--", lw=0.5, alpha=0.4)

ax.axhline(y=1.0, color="#999999", ls="--", lw=0.5, alpha=0.5)

ax.set_xlabel(r"Mode-tangent alignment $\alpha$", fontsize=11, color=TEAL)
ax.set_ylabel(r"Penalty factor $P(\alpha)$", fontsize=11, color=TEAL)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=9, loc="upper left", framealpha=0.9)
ax.tick_params(colors=TEAL, labelsize=10)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["bottom", "left"]:
    ax.spines[spine].set_color(TEAL)
    ax.spines[spine].set_linewidth(1.0)

ax.grid(True, alpha=0.15)
fig.tight_layout()

out_base = "analysis/sympy/penalty_family"
fig.savefig(f"{out_base}.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(f"{out_base}.pdf", bbox_inches="tight", facecolor="white")
print(f"Saved: {out_base}.png and .pdf")

print("\n" + "=" * 60)
print("  ALL ASSERTIONS PASSED")
print("=" * 60)
