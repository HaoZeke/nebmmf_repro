Require Import Reals.
Require Import Lra.
Open Scope R_scope.

(* b = 1/(1+s) is in (0,1) for s > 0 *)
Theorem derived_base_in_unit_interval :
  forall s : R, 0 < s ->
    let b := 1 / (1 + s) in
    0 < b /\ b < 1.
Proof.
  intros s Hs. split.
  - unfold Rdiv. apply Rmult_lt_0_compat. lra. apply Rinv_0_lt_compat. lra.
  - apply Rmult_lt_reg_r with (r := 1 + s). lra.
    unfold Rdiv. rewrite Rmult_1_l.
    rewrite Rinv_l. lra. lra.
Qed.

(* At the default strength s=1.5, the derived base is 0.4 *)
(* This is 1/(1+1.5) = 1/2.5 = 2/5 *)
Theorem default_base_value :
  1 / (1 + 3/2) = 2/5.
Proof.
  field.
Qed.
