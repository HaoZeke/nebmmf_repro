Require Import Reals.
Require Import Lra.
Open Scope R_scope.

(* If T_adaptive >= T_absolute, then
   (F < T_adaptive \/ F < T_absolute) <-> F < T_adaptive *)
Theorem trigger_force_redundant :
  forall F T_adaptive T_absolute : R,
    T_adaptive >= T_absolute ->
    (F < T_adaptive \/ F < T_absolute) <-> F < T_adaptive.
Proof.
  intros F Ta Tf Hge. split.
  - intros [H | H]. exact H. lra.
  - intros H. left. exact H.
Qed.

(* Under default parameters: trigger_force = 2*f_tol, so
   T_adaptive >= 2*f_tol implies T_adaptive >= trigger_force *)
Theorem default_params_redundant :
  forall f_tol T_adaptive : R,
    0 < f_tol ->
    T_adaptive >= 2 * f_tol ->
    let trigger_force := 2 * f_tol in
    T_adaptive >= trigger_force.
Proof.
  intros. unfold trigger_force. lra.
Qed.
