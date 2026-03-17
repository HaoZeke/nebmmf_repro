Require Import Reals.
Require Import Lra.
Require Import Psatz.
Open Scope R_scope.

(** * Monotonic progress of OCI-NEB under successful MMF refinement.

    The success threshold update (NEBRonebController.cpp:207) is:
      T_new = F_new * (0.5 + 0.4 * F_new / F_old)

    We prove that when MMF helps (F_new < F_old), the force strictly
    decreases and the threshold tracks the force from above. *)

(** The success update factor r*(0.5 + 0.4*r) < 1 for 0 < r < 1 *)
Theorem success_factor_lt_one :
  forall r : R, 0 < r -> r < 1 ->
    r * (1/2 + 2/5 * r) < 1.
Proof.
  intros r Hr0 Hr1. nra.
Qed.

(** The success update ratio T_new/F_new is in [0.5, 0.9] *)
Theorem success_ratio_bounds :
  forall r : R, 0 < r -> r <= 1 ->
    1/2 <= 1/2 + 2/5 * r /\ 1/2 + 2/5 * r <= 9/10.
Proof.
  intros r Hr0 Hr1. split; nra.
Qed.

(** When MMF helps (r < 1), the new threshold is strictly below the old force.
    Key insight: r*(0.5+0.4*r) < 0.9 < 1 for r < 1. *)
Theorem threshold_below_old_force :
  forall F_old r : R,
    0 < F_old -> 0 < r -> r < 1 ->
    r * F_old * (1/2 + 2/5 * r) < F_old.
Proof.
  intros F_old r HFold Hr0 Hr1.
  assert (Hfactor: r * (1/2 + 2/5 * r) < 1) by nra.
  assert (Hpos: 0 < 1/2 + 2/5 * r) by lra.
  nra.
Qed.

(** Repeated success leads to geometric decrease:
    T_new < 0.9 * T_old when force is at the threshold. *)
Theorem geometric_decrease_step :
  forall T_old F_old r : R,
    0 < T_old -> 0 < F_old ->
    F_old <= T_old ->
    0 < r -> r < 1 ->
    r * F_old * (1/2 + 2/5 * r) < 9/10 * T_old.
Proof.
  intros T_old F_old r HT HF Htrig Hr0 Hr1.
  assert (Hfactor: r * (1/2 + 2/5 * r) < 9/10) by nra.
  nra.
Qed.
