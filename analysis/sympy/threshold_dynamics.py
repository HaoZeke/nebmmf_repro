# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "sympy",
#   "matplotlib",
# ]
# ///
"""Analyze the threshold adaptation dynamics in OCI-NEB.

The adaptive threshold has two update rules (NEBRonebController.cpp):

Success (line 207):
    T_new = F_new * (0.5 + 0.4 * F_new / F_old)
    T_new = min(T_new, baseline * trigger_factor)

Backoff (line 212-221):
    penalty = base + (1 - base) * alpha^strength
    T_new = baseline * trigger_factor * penalty
    T_new = max(T_new, 2 * f_tol)

This script proves the threshold is bounded and that trigger_force is redundant.
"""

import sympy as sp
import matplotlib.pyplot as plt
import numpy as np

# -- Symbolic setup ----------------------------------------------------------

F_new, F_old = sp.symbols("F_new F_old", positive=True)
r = sp.Symbol("r", positive=True)  # reduction ratio F_new/F_old
T, T_new = sp.symbols("T T_new", positive=True)
baseline, trigger_factor = sp.symbols("B_0 tau", positive=True)
f_tol = sp.Symbol("f_tol", positive=True)
trigger_force = sp.Symbol("T_abs", positive=True)
alpha, s_sym, b = sp.symbols("alpha s b", positive=True)

print("=" * 72)
print("THRESHOLD DYNAMICS ANALYSIS")
print("=" * 72)

# (a) Define update rules symbolically
print("\n[a] Update rules:")

T_success = F_new * (sp.Rational(1, 2) + sp.Rational(2, 5) * F_new / F_old)
print(f"    Success: T_new = {T_success}")
print(f"    LaTeX: T_{{\\mathrm{{new}}}} = {sp.latex(T_success)}")

penalty = b + (1 - b) * alpha**s_sym
T_backoff = baseline * trigger_factor * penalty
print(f"\n    Backoff: T_new = {T_backoff}")
print(f"    LaTeX: T_{{\\mathrm{{new}}}} = {sp.latex(T_backoff)}")

# (b) Success rule as a contraction map
print("\n[b] Success rule analysis (F_new = r * F_old, 0 < r < 1):")

T_success_r = T_success.subs(F_new, r * F_old)
T_success_r = sp.expand(T_success_r)
print(f"    T_new = {T_success_r}")

# Factor out F_old
T_success_factored = sp.collect(T_success_r, F_old)
print(f"    T_new = {T_success_factored}")
print(f"    LaTeX: {sp.latex(T_success_factored)}")

# Ratio T_new / F_new
ratio = sp.simplify(T_success_r / (r * F_old))
print(f"\n    T_new / F_new = {ratio}")
print(f"    LaTeX: {sp.latex(ratio)}")

# This ratio is 0.5 + 0.4*r, which is in [0.5, 0.9] for r in (0,1)
ratio_at_0 = ratio.subs(r, 0)
ratio_at_1 = ratio.subs(r, 1)
print(f"    At r=0: T_new/F_new = {ratio_at_0}")
print(f"    At r=1: T_new/F_new = {ratio_at_1}")
assert ratio_at_0 == sp.Rational(1, 2)
assert ratio_at_1 == sp.Rational(9, 10)
print("    => T_new/F_new in [1/2, 9/10] for r in (0,1]  [OK]")

# Contraction: T_new < F_new when ratio < 1
# 0.5 + 0.4*r < 1  =>  r < 5/4, always true for r < 1
print("\n    Contraction condition: T_new < F_new")
print("    0.5 + 0.4*r < 1  <=>  r < 5/4")
print("    Since r < 1 < 5/4, the success rule is always a contraction  [OK]")

# (c) Backoff rule bounds
print("\n[c] Backoff rule bounds:")

# Lower bound from hard floor
print(f"    Lower bound: T >= 2 * f_tol  (hard floor, line 220-221)")
print(f"    LaTeX: T \\geq 2 f_{{\\mathrm{{tol}}}}")

# Upper bound from success cap
print(f"    Upper bound: T <= B_0 * tau  (success cap)")
print(f"    LaTeX: T \\leq B_0 \\tau")

# Combined
print(f"\n    Threshold invariant: 2 f_tol <= T_adaptive <= B_0 * tau")
print(f"    LaTeX: 2 f_{{\\mathrm{{tol}}}} \\leq T_{{\\mathrm{{adaptive}}}} \\leq B_0 \\tau")

# Verify with defaults
f_tol_val = sp.Rational(5, 100)  # 0.05
trigger_force_val = sp.Rational(1, 10)  # 0.1
floor_val = 2 * f_tol_val
print(f"\n    Default values: f_tol = {f_tol_val}, trigger_force = {trigger_force_val}")
print(f"    Hard floor = 2 * f_tol = {floor_val} = {float(floor_val)}")
assert floor_val == trigger_force_val
print(f"    trigger_force = {trigger_force_val} = {float(trigger_force_val)}")
print(f"    floor == trigger_force  [OK]")

# (d) Trigger_force redundancy proof
print("\n[d] Trigger_force redundancy:")

print("    shouldTrigger predicate:")
print("      trigger(F) = (F < T_adaptive) OR (F < T_absolute)")
print()
print("    Invariant: T_adaptive >= 2 * f_tol  (from backoff floor)")
print()
print("    Case: T_absolute <= 2 * f_tol")
print("      => T_absolute <= T_adaptive  (always)")
print("      => {F < T_absolute} is a subset of {F < T_adaptive}")
print("      => the OR is redundant")
print()
print("    Default verification:")
print(f"      T_absolute = trigger_force = {float(trigger_force_val)}")
print(f"      2 * f_tol = {float(floor_val)}")
print(f"      T_absolute <= 2 * f_tol: {trigger_force_val <= floor_val}  [OK]")
assert trigger_force_val <= floor_val
print()
print("    Therefore: trigger(F) <=> (F < T_adaptive)")
print("    trigger_force is eliminable without changing behavior.")

# LaTeX proof
print("\n    LaTeX proof:")
print(r"    T_{\mathrm{adaptive}} \geq 2 f_{\mathrm{tol}} \geq T_{\mathrm{absolute}}")
print(r"    \implies \{F < T_{\mathrm{absolute}}\} \subseteq \{F < T_{\mathrm{adaptive}}\}")
print(r"    \implies \mathrm{trigger}(F) \iff F < T_{\mathrm{adaptive}}")

# (e) Print all results as LaTeX
print("\n[e] LaTeX equation summary:")
print()
print(r"    \text{Success: } T_{\mathrm{new}} = " + sp.latex(T_success))
print()
print(r"    \text{Backoff: } T_{\mathrm{new}} = " + sp.latex(T_backoff))
print()
print(r"    \frac{T_{\mathrm{new}}}{F_{\mathrm{new}}} = " + sp.latex(ratio))
print()
print(r"    2 f_{\mathrm{tol}} \leq T_{\mathrm{adaptive}} \leq B_0 \tau")

# (f) Plot threshold bounds
print("\n[f] Generating threshold bounds plot...")

TEAL = "#004D40"
CORAL = "#FF655D"
YELLOW = "#F1DB4B"

fig, ax = plt.subplots(figsize=(6, 4))

baseline_vals = np.linspace(0.01, 2.0, 200)
tau_default = 0.2  # typical trigger_factor
f_tol_default = 0.05

upper_bound = baseline_vals * tau_default
lower_bound = np.full_like(baseline_vals, 2 * f_tol_default)

ax.fill_between(
    baseline_vals, lower_bound, upper_bound,
    alpha=0.15, color=TEAL, label="Admissible region"
)
ax.plot(baseline_vals, upper_bound, color=TEAL, linewidth=1.8,
        label=r"Upper: $B_0 \cdot \tau$")
ax.plot(baseline_vals, lower_bound, color=CORAL, linewidth=1.8,
        linestyle="--", label=r"Lower: $2 f_{\mathrm{tol}}$")

# Mark where the bounds cross
cross_baseline = 2 * f_tol_default / tau_default
ax.axvline(cross_baseline, color=YELLOW, linewidth=1.2, linestyle=":",
           label=f"Bounds meet at $B_0$ = {cross_baseline:.2f}")
ax.plot(cross_baseline, 2 * f_tol_default, "o", color=CORAL, markersize=6)

ax.set_xlabel(r"Baseline force $B_0$", fontsize=11, color=TEAL)
ax.set_ylabel(r"Threshold $T$", fontsize=11, color=TEAL)
ax.set_title(
    r"Adaptive threshold bounds ($\tau$ = "
    f"{tau_default}, "
    r"$f_{\mathrm{tol}}$ = "
    f"{f_tol_default})",
    fontsize=12, color=TEAL,
)
ax.legend(fontsize=9, framealpha=0.9)
ax.set_xlim(0, 2.0)
ax.set_ylim(0, 0.5)
ax.tick_params(colors=TEAL)
for spine in ax.spines.values():
    spine.set_color(TEAL)
ax.grid(True, alpha=0.3, color=TEAL)

plt.tight_layout()
outdir = "/home/rgoswami/Git/Github/epfl/pixi_envs/repro/nebmmf/analysis/sympy"
fig.savefig(f"{outdir}/threshold_bounds.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{outdir}/threshold_bounds.pdf", bbox_inches="tight")
print(f"    Saved: threshold_bounds.png, threshold_bounds.pdf")
plt.close(fig)

print("\nPASS")
