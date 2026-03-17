Require Import Reals.
Require Import Lra.
Open Scope R_scope.

(* The backoff update always produces T >= 2*f_tol by construction *)
Theorem threshold_lower_bound :
  forall (f_tol raw_threshold : R),
    0 < f_tol ->
    Rmax raw_threshold (2 * f_tol) >= 2 * f_tol.
Proof.
  intros. unfold Rmax. destruct (Rle_dec raw_threshold (2 * f_tol)); lra.
Qed.
