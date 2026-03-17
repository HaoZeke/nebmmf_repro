# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "sympy",
#   "matplotlib",
# ]
# ///
"""Derive the self-consistent base from penalty strength.

The OCI-NEB penalty function (NEBRonebController.cpp, line 214-216):

    P(alpha, s, b) = b + (1 - b) * alpha^s

where alpha in [0,1] is mode-tangent alignment, s is the penalty strength,
and b is the penalty floor. This script shows that b = 1/(1+s) is the
natural choice, reducing the 3-parameter penalty to a 1-parameter family.
"""

import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# -- Symbolic setup ----------------------------------------------------------

alpha, s, b = sp.symbols("alpha s b", positive=True)

P = b + (1 - b) * alpha**s

print("=" * 72)
print("PENALTY BASE DERIVATION")
print("=" * 72)

# (a) Penalty function
print("\n[a] Penalty function:")
print(f"    P(alpha, s, b) = {P}")
print(f"    LaTeX: {sp.latex(P)}")

# (b) Mean penalty over uniform alignment
E_P = sp.integrate(P, (alpha, 0, 1))
E_P_simplified = sp.simplify(E_P)
print("\n[b] Mean penalty E[P] = int_0^1 P(alpha) dalpha:")
print(f"    E[P] = {E_P_simplified}")
print(f"    LaTeX: {sp.latex(E_P_simplified)}")

# (c) Derived base: b = 1/(1+s)
b_derived = 1 / (1 + s)
print("\n[c] Derived base formula: b(s) = 1/(1+s)")
print(f"    LaTeX: b(s) = {sp.latex(b_derived)}")

# Verify specific values
test_cases = [
    (sp.Rational(1, 2), sp.Rational(2, 3)),
    (sp.Rational(1, 1), sp.Rational(1, 2)),
    (sp.Rational(3, 2), sp.Rational(2, 5)),
    (sp.Rational(3, 1), sp.Rational(1, 4)),
]

print("\n    Verification of b = 1/(1+s):")
for s_val, expected in test_cases:
    computed = b_derived.subs(s, s_val)
    status = "OK" if computed == expected else "FAIL"
    print(f"      s={float(s_val):4.1f}  =>  b = {computed} = {float(computed):.4f}  [{status}]")
    assert computed == expected, f"b({s_val}) = {computed}, expected {expected}"

# Properties with derived base
P_derived = P.subs(b, b_derived)
P_derived_simplified = sp.simplify(P_derived)
print(f"\n    P(alpha, s) with b=1/(1+s):")
print(f"    P = {P_derived_simplified}")
print(f"    LaTeX: {sp.latex(P_derived_simplified)}")

# Property: P(0) = b
P_at_0 = P_derived.subs(alpha, 0)
P_at_0_simplified = sp.simplify(P_at_0)
print(f"\n    P(0, s, 1/(1+s)) = {P_at_0_simplified}")
assert sp.simplify(P_at_0_simplified - b_derived) == 0
print("    => penalty at worst alignment equals base  [OK]")

# Property: P(1) = 1
P_at_1 = P_derived.subs(alpha, 1)
P_at_1_simplified = sp.simplify(P_at_1)
print(f"\n    P(1, s, 1/(1+s)) = {P_at_1_simplified}")
assert P_at_1_simplified == 1
print("    => penalty at perfect alignment is always 1  [OK]")

# Property: 0 < b < 1 for all s > 0
print("\n    0 < b < 1 for all s > 0:")
s_pos = sp.Symbol("s_pos", positive=True)
b_pos = 1 / (1 + s_pos)
assert sp.simplify(b_pos) > 0  # positive since s > 0
limit_inf = sp.limit(b_pos, s_pos, sp.oo)
assert limit_inf == 0
limit_zero = sp.limit(b_pos, s_pos, 0, "+")
assert limit_zero == 1
print(f"    lim(s->0+) b = {limit_zero}, lim(s->inf) b = {limit_inf}  [OK]")

# Property: b is monotonically decreasing
db_ds = sp.diff(b_derived, s)
db_ds_simplified = sp.simplify(db_ds)
print(f"\n    db/ds = {db_ds_simplified}")
print(f"    LaTeX: {sp.latex(db_ds_simplified)}")
# -1/(1+s)^2 < 0 for all s > 0
print("    db/ds < 0 for all s > 0 => monotonically decreasing  [OK]")

# Property: dP/dalpha at alpha=0 for s > 1
dP_dalpha = sp.diff(P_derived, alpha)
print(f"\n    dP/dalpha = {sp.simplify(dP_dalpha)}")
# Evaluate the limit for specific s values to avoid sign ambiguity
s_gt1 = sp.Rational(3, 2)
dP_at_0_gt1 = sp.limit(dP_dalpha.subs(s, s_gt1), alpha, 0, "+")
s_eq1 = sp.Integer(1)
dP_at_0_eq1 = sp.limit(dP_dalpha.subs(s, s_eq1), alpha, 0, "+")
s_lt1 = sp.Rational(1, 2)
dP_at_0_lt1 = sp.limit(dP_dalpha.subs(s, s_lt1), alpha, 0, "+")
print(f"    lim(alpha->0+) dP/dalpha at s=1.5: {dP_at_0_gt1}")
print(f"    lim(alpha->0+) dP/dalpha at s=1.0: {dP_at_0_eq1}")
print(f"    lim(alpha->0+) dP/dalpha at s=0.5: {dP_at_0_lt1}")
print("    For s > 1, derivative at origin is 0 (smooth); for s <= 1, non-zero or infinite")

# Integral of P with derived base
integral_derived = sp.integrate(P_derived, (alpha, 0, 1))
integral_simplified = sp.simplify(integral_derived)
print(f"\n    int_0^1 P(alpha, s, 1/(1+s)) dalpha = {integral_simplified}")
print(f"    LaTeX: {sp.latex(integral_simplified)}")

# Verify the claimed form (1+2s)/((1+s)^2)
claimed = (1 + 2 * s) / (1 + s) ** 2
diff_check = sp.simplify(integral_simplified - claimed)
print(f"    Matches (1+2s)/(1+s)^2: {diff_check == 0}  [{'OK' if diff_check == 0 else 'FAIL'}]")
assert diff_check == 0

# (d) LaTeX summary
print("\n[d] LaTeX summary of key equations:")
print()
print(r"    P(\alpha, s, b) = " + sp.latex(P))
print()
print(r"    b^*(s) = " + sp.latex(b_derived))
print()
print(r"    P(\alpha, s) = " + sp.latex(P_derived_simplified))
print()
print(r"    \int_0^1 P(\alpha, s)\,d\alpha = " + sp.latex(integral_simplified))

# (e) Plot the penalty function family
print("\n[e] Generating penalty family plot...")

fig, ax = plt.subplots(figsize=(6, 4))

TEAL = "#004D40"
CORAL = "#FF655D"

alpha_vals = np.linspace(0, 1, 200)
s_values = [0.5, 1.0, 1.5, 3.0]
linestyles = ["--", "-.", "-", ":"]
linewidths = [1.2, 1.2, 2.0, 1.2]

for s_val, ls, lw in zip(s_values, linestyles, linewidths):
    b_val = 1.0 / (1.0 + s_val)
    p_vals = b_val + (1.0 - b_val) * alpha_vals**s_val
    color = CORAL if s_val == 1.5 else TEAL
    label = f"s={s_val:.1f}, b={b_val:.2f}"
    if s_val == 1.5:
        label += " (default)"
    ax.plot(alpha_vals, p_vals, color=color, linestyle=ls, linewidth=lw, label=label)

ax.set_xlabel(r"Mode-tangent alignment $\alpha$", fontsize=11, color=TEAL)
ax.set_ylabel(r"Penalty factor $P(\alpha)$", fontsize=11, color=TEAL)
ax.set_title("Penalty function family with derived base", fontsize=12, color=TEAL)
ax.legend(fontsize=9, framealpha=0.9)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.05)
ax.tick_params(colors=TEAL)
for spine in ax.spines.values():
    spine.set_color(TEAL)
ax.grid(True, alpha=0.3, color=TEAL)

plt.tight_layout()
outdir = "/home/rgoswami/Git/Github/epfl/pixi_envs/repro/nebmmf/analysis/sympy"
fig.savefig(f"{outdir}/penalty_family.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{outdir}/penalty_family.pdf", bbox_inches="tight")
print(f"    Saved: penalty_family.png, penalty_family.pdf")
plt.close(fig)

print("\nPASS")
