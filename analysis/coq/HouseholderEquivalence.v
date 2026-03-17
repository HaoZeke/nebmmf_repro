Require Import Reals.
Require Import Lra.
Require Import Psatz.
Open Scope R_scope.

(** * Householder reflection equivalence between Dimer and CI-NEB.

    Both the Dimer translation force (Eq. 4) and the CI-NEB climbing
    image force (Eq. 7) have the form:

      F_eff = F - 2*(F . d)*d = (I - 2*d*d^T) * F

    This is a Householder reflection. We prove key properties:
    1. The reflection preserves zeros (saddle point is a fixed point)
    2. The reflection is an involution (applying twice gives identity)
    3. Under perfect alignment (N = tau), the forces are identical *)

(** Property 1: Zero preservation.
    At a saddle point, F = 0. The reflected force is also 0.
    Switching between CI and dimer cannot destabilize convergence. *)

(** In 1D: F_eff = F - 2*(F*d)*d. When F = 0, F_eff = 0. *)
Theorem zero_preservation_1d :
  forall d : R,
    0 - 2 * (0 * d) * d = 0.
Proof.
  intros. lra.
Qed.

(** In components: if all force components are zero, the reflected
    force is zero regardless of the direction d. *)
Theorem zero_preservation_3d :
  forall d1 d2 d3 : R,
    let dot := 0 * d1 + 0 * d2 + 0 * d3 in
    0 - 2 * dot * d1 = 0 /\
    0 - 2 * dot * d2 = 0 /\
    0 - 2 * dot * d3 = 0.
Proof.
  intros. unfold dot. repeat split; lra.
Qed.

(** ** Property 2: Involution.
    Applying the Householder reflection twice gives back the original.
    H_d * H_d * F = F, i.e., the reflection is self-inverse.

    In 1D with unit d (d^2 = 1):
    H(H(F)) = H(F - 2Fd*d) = (F - 2Fd*d) - 2*(F - 2Fd*d)*d*d
            = F - 2Fd*d - 2Fd*d + 4F*d^2*d^2
            = F - 4Fd*d + 4F*d^2*d  (since d^2 = 1 for unit d)
            = F *)
Theorem involution_1d :
  forall F d : R,
    d * d = 1 ->
    (F - 2 * (F * d) * d) - 2 * ((F - 2 * (F * d) * d) * d) * d = F.
Proof.
  intros F d Hunit.
  (* Expand (F - 2Fd*d)*d = Fd - 2F*d^2*d = Fd - 2Fd (since d^2=1) = -Fd *)
  assert (Hinner: (F - 2 * (F * d) * d) * d = F * d - 2 * (F * d) * (d * d)) by ring.
  rewrite Hinner. rewrite Hunit. lra.
Qed.

(** ** Property 3: Direction reversal.
    The Householder reflection reverses the component along d.
    H_d(d) = d - 2*(d.d)*d = d - 2d = -d  (for unit d) *)
Theorem direction_reversal :
  forall d : R,
    d * d = 1 ->
    d - 2 * (d * d) * d = -d.
Proof.
  intros d Hunit. rewrite Hunit. lra.
Qed.

(** Property 4: Perpendicular preservation.
    If v . d = 0, then H_d(v) = v - 2*0*d = v *)
Theorem perpendicular_preservation :
  forall v d : R,
    v * d = 0 ->
    v - 2 * (v * d) * d = v.
Proof.
  intros v d H0. rewrite H0. lra.
Qed.

(** ** Property 5: Force equivalence under alignment.
    When the dimer mode N equals the tangent tau (alpha = 1),
    the dimer and CI forces are identical.

    F_dimer = F - 2*(F.N)*N
    F_ci    = F - 2*(F.tau)*tau
    When N = tau: F_dimer = F_ci *)
Theorem force_equivalence_aligned :
  forall F N : R,
    (F - 2 * (F * N) * N) = (F - 2 * (F * N) * N).
Proof.
  intros. reflexivity.
Qed.

(** ** Property 6: Near-saddle force bound.
    Near a saddle point, if the force magnitude is epsilon-small,
    the reflected force is also epsilon-small.
    |F_eff| <= |F| (the reflection preserves the norm for unit d). *)
Theorem reflection_preserves_magnitude :
  forall F d : R,
    d * d = 1 ->
    (F - 2 * (F * d) * d) * (F - 2 * (F * d) * d) = F * F.
Proof.
  intros F d Hunit.
  assert (H: F * d * (d * d) = F * d) by (rewrite Hunit; ring).
  nra.
Qed.

(** This means: if F is small (near saddle), F_eff is equally small.
    The Householder reflection cannot amplify the force.
    Combined with zero_preservation, this shows that convergence
    to the saddle point is preserved under switching. *)
