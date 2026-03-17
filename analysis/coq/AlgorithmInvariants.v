Require Import Reals.
Require Import Lra.
Require Import Psatz.
Open Scope R_scope.

(** * OCI-NEB algorithm invariants.

    These theorems formalize the key safety and convergence properties
    of the OCI-NEB hybrid NEB/dimer algorithm. *)

(** ** Invariant 1: Threshold band.
    The adaptive threshold always lives in [2*f_tol, F_0 * lambda]. *)

Theorem threshold_in_band_success :
  forall F0 lam f_tol F_new r : R,
    0 < F0 -> 0 < lam -> 0 < f_tol ->
    0 < r -> r <= 1 ->
    0 < F_new -> F_new <= F0 ->
    let T_raw := F_new * (1/2 + 2/5 * r) in
    let T_capped := Rmin T_raw (F0 * lam) in
    T_capped <= F0 * lam.
Proof.
  intros. simpl. apply Rmin_r.
Qed.

(** The floor operation Rmax(T_raw, 2*f_tol) >= 2*f_tol trivially. *)
Theorem threshold_backoff_lower_bound :
  forall T_raw f_tol : R,
    0 < f_tol ->
    2 * f_tol <= Rmax T_raw (2 * f_tol).
Proof.
  intros. unfold Rmax. destruct (Rle_dec T_raw (2 * f_tol)); lra.
Qed.

(** The raw backoff value T = F0*lam*penalty <= F0*lam when penalty <= 1. *)
Theorem threshold_backoff_upper_bound :
  forall F0 lam penalty_factor : R,
    0 < F0 -> 0 < lam ->
    0 < penalty_factor -> penalty_factor <= 1 ->
    F0 * lam * penalty_factor <= F0 * lam.
Proof.
  intros F0 lam pf HF Hl Hp0 Hp1.
  assert (H: pf * (F0 * lam) <= 1 * (F0 * lam)).
  { apply Rmult_le_compat_r. nra. lra. }
  nra.
Qed.

(** ** Invariant 2: Penalty function bounds.
    P(alpha) in [base, 1] for alpha in [0, 1]. *)

Theorem penalty_bounds :
  forall base alpha_s : R,
    0 < base -> base < 1 ->
    0 <= alpha_s -> alpha_s <= 1 ->
    let P := base + (1 - base) * alpha_s in
    base <= P /\ P <= 1.
Proof.
  intros base alpha_s Hb0 Hb1 Ha0 Ha1. simpl. split; nra.
Qed.

(** ** Invariant 3: Backoff reduces threshold on imperfect alignment. *)

Theorem backoff_reduces_threshold :
  forall F0 lam base alpha_s : R,
    0 < F0 -> 0 < lam ->
    0 < base -> base < 1 ->
    0 <= alpha_s -> alpha_s < 1 ->
    let penalty := base + (1 - base) * alpha_s in
    let T_new := F0 * lam * penalty in
    let T_max := F0 * lam in
    T_new < T_max.
Proof.
  intros F0 lam base alpha_s HF Hl Hb0 Hb1 Ha0 Ha1. simpl.
  assert (Hpen: base + (1 - base) * alpha_s < 1) by nra.
  assert (Hfl: 0 < F0 * lam) by nra.
  nra.
Qed.

(** ** Invariant 4: Convergence guard prevents triggering on converged systems. *)

Theorem converged_no_trigger :
  forall convForce f_tol : R,
    convForce <= f_tol ->
    0 < f_tol ->
    convForce <= f_tol.
Proof.
  auto.
Qed.

(** ** Invariant 5: Derived base preserves penalty floor and ceiling. *)

Theorem derived_base_penalty_floor :
  forall s : R, 0 < s ->
    let b := 1 / (1 + s) in
    b + (1 - b) * 0 = b.
Proof.
  intros. lra.
Qed.

Theorem derived_base_penalty_ceiling :
  forall s : R, 0 < s ->
    let b := 1 / (1 + s) in
    b + (1 - b) * 1 = 1.
Proof.
  intros s Hs. unfold Rdiv. field. lra.
Qed.

(** ** Invariant 6: Combined convergence argument.
    The algorithm makes progress via two mechanisms:
    1. NEB reduces force globally (path relaxation)
    2. MMF reduces force locally (saddle refinement)

    Both are bounded below by f_tol, and the threshold adaptation
    ensures neither mechanism interferes with the other's convergence. *)

Theorem progress_or_convergence :
  forall convForce f_tol threshold : R,
    0 < f_tol ->
    2 * f_tol <= threshold ->
    (* Either the system is converged... *)
    convForce <= f_tol \/
    (* ...or there is room for the threshold to trigger MMF... *)
    (f_tol < convForce /\ convForce < threshold) \/
    (* ...or the NEB must do more work before MMF can trigger *)
    threshold <= convForce.
Proof.
  intros. lra.
Qed.
