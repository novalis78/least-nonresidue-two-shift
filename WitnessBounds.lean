import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity

/-! Analytic tail estimate and final comparison for the two-shift certificate.
The interval enclosures on the finite domains are supplied by the Python verifiers. -/
noncomputable section
open MeasureTheory Set
namespace TwoShift

/-- A summand of the spatial amplitude. Its denominator has degree m+2. -/
def term (b ω : ℝ) (m : ℕ) (z : ℝ) : ℂ :=
  (b : ℂ) * Complex.exp (((ω*z : ℝ) : ℂ)*Complex.I) /
    (1 + (z : ℂ)*Complex.I)^(m+2)

theorem term_bound (b ω : ℝ) (m : ℕ) {z : ℝ} (hz : 0 < z) :
    ‖term b ω m z‖ ≤ |b| / z^(m+2) := by
  have h : z ≤ ‖1 + (z : ℂ)*Complex.I‖ := by
    simpa using Complex.im_le_norm (1 + (z : ℂ)*Complex.I)
  rw [term, norm_div, norm_mul, Complex.norm_exp_ofReal_mul_I, mul_one, norm_pow]
  rw [Complex.norm_real, Real.norm_eq_abs]
  exact div_le_div_of_nonneg_left (abs_nonneg _) (pow_pos hz _) (pow_le_pow_left₀ hz.le h _)

private theorem envelope_eq (b : ℝ) (m : ℕ) {x : ℝ} (hx : 0 < x) :
    |b| / x^(m+2) = |b| * x^(-((m+2 : ℕ) : ℝ)) := by
  rw [Real.rpow_neg hx.le, Real.rpow_natCast, div_eq_mul_inv]

private theorem envelope_integrable (b : ℝ) (m : ℕ) {R : ℝ} (hR : 0 < R) :
    IntegrableOn (fun z : ℝ => |b| / z^(m+2)) (Ioi R) := by
  have hk : -((m+2 : ℕ) : ℝ) < -1 := by push_cast; linarith
  exact IntegrableOn.congr_fun ((integrableOn_Ioi_rpow_of_lt hk hR).const_mul |b|)
    (fun x hx => (envelope_eq b m (hR.trans hx)).symm) measurableSet_Ioi

theorem envelope_integral (b : ℝ) (m : ℕ) {R : ℝ} (hR : 0 < R) :
    (∫ z : ℝ in Ioi R, |b| / z^(m+2)) = |b| / ((m+1 : ℝ)*R^(m+1)) := by
  have hk : -((m+2 : ℕ) : ℝ) < -1 := by push_cast; linarith
  calc
    _ = ∫ z : ℝ in Ioi R, |b| * z^(-((m+2 : ℕ) : ℝ)) :=
      setIntegral_congr_fun measurableSet_Ioi (fun x hx => envelope_eq b m (hR.trans hx))
    _ = _ := by
      rw [integral_const_mul, integral_Ioi_rpow_of_lt hk hR]
      rw [show -((m+2 : ℕ) : ℝ) + 1 = -((m+1 : ℕ) : ℝ) by push_cast; ring]
      rw [Real.rpow_neg hR.le, Real.rpow_natCast]
      push_cast
      rw [neg_div_neg_eq]
      simp [div_eq_mul_inv, mul_comm]

/-- The contribution beyond R is integrable and bounded by the finite
    sum used in the verifiers, before the common factor 1/π. -/
theorem tail_bound {ι : Type*} (s : Finset ι) (b ω : ι → ℝ) (m : ι → ℕ)
    {R : ℝ} (hR : 0 < R) :
    IntegrableOn (fun z => ‖∑ j ∈ s, term (b j) (ω j) (m j) z‖) (Ioi R) ∧
    (∫ z : ℝ in Ioi R, ‖∑ j ∈ s, term (b j) (ω j) (m j) z‖) ≤
      ∑ j ∈ s, |b j| / ((m j+1 : ℝ)*R^(m j+1)) := by
  have env := fun j => envelope_integrable (b j) (m j) hR
  have sumInt := integrable_finset_sum s (fun j _ => env j)
  have pointwise (z : ℝ) (hz : z ∈ Ioi R) :
      ‖∑ j ∈ s, term (b j) (ω j) (m j) z‖ ≤ ∑ j ∈ s, |b j| / z^(m j+2) :=
    (norm_sum_le _ _).trans (Finset.sum_le_sum (fun j _ => term_bound _ _ _ (hR.trans hz)))
  have meas : Measurable (fun z => ‖∑ j ∈ s, term (b j) (ω j) (m j) z‖) := by
    unfold term
    fun_prop
  have hi : IntegrableOn (fun z => ‖∑ j ∈ s, term (b j) (ω j) (m j) z‖) (Ioi R) := by
    apply sumInt.mono' meas.aestronglyMeasurable
    filter_upwards [ae_restrict_mem (μ := volume) measurableSet_Ioi] with z hz
    simpa using pointwise z hz
  refine ⟨hi, ?_⟩
  calc
    _ ≤ ∫ z : ℝ in Ioi R, ∑ j ∈ s, |b j| / z^(m j+2) :=
      setIntegral_mono_on hi sumInt measurableSet_Ioi pointwise
    _ = _ := by
      rw [integral_finset_sum s (fun j _ => env j)]
      exact Finset.sum_congr rfl (fun j _ => envelope_integral (b j) (m j) hR)

/-- Integer numerators from Table 1, read by rows of polynomial degree. -/
def coefficients : List (ℤ × ℤ) :=
  [(1295, 0), (79929, -9074),
   (-6609595, 5043626), (117992958, -14383101),
   (-234402020, -118792841), (332397224, 61384470),
   (-46679113, -118505260), (20002165, 154829818),
   (-2996805, -138162493), (-22963877, 320107061),
   (279853645, 234086697), (-809901785, 62986550)]

/-- The exact rational right-hand side of the tail estimate at R = 64. -/
def witnessTail : ℚ :=
  (coefficients.zipIdx.map fun ((b₁, b₂), j) =>
    ((b₁.natAbs + b₂.natAbs : ℕ) : ℚ) /
      (10^7 * (2*(j:ℚ)+1) * 64^(2*j+1))).sum

theorem witness_tail_arithmetic : witnessTail < 2035/1000000000 := by decide +kernel

/-- The two rational bounds are rounded outward from the Arb output. -/
theorem ratio_of_enclosures {D N : ℝ}
    (hD₀ : 0 ≤ D) (hD : D ≤ 0.8005327903879) (hN : 0.9177970191817 ≤ N) :
    (D/N)^2 < 0.76079030 := by
  have hN₀ : 0 < N := lt_of_lt_of_le (by norm_num) hN
  have hr : D/N ≤ 0.8005327903879 / 0.9177970191817 :=
    div_le_div₀ (by norm_num) hD (by norm_num) hN
  exact lt_of_le_of_lt (pow_le_pow_left₀ (div_nonneg hD₀ hN₀.le) hr 2) (by norm_num)

/-- Transfer the real enclosures to any positive-score upper bound C.
For the Fourier extremal constant the last hypothesis is the supremum property. -/
theorem constant_of_enclosures {C D N : ℝ}
    (hD₀ : 0 < D) (hD : D ≤ 0.8005327903879) (hN : 0.9177970191817 ≤ N)
    (hscore : N/D ≤ C) : (1/C)^2 < 0.76079030 := by
  have hN₀ : 0 < N := lt_of_lt_of_le (by norm_num) hN
  have hC : 0 < C := (div_pos hN₀ hD₀).trans_le hscore
  have hInv : 1/C ≤ D/N := by
    apply (div_le_div_iff₀ hC hN₀).2
    simpa [mul_comm] using (div_le_iff₀ hD₀).1 hscore
  exact (pow_le_pow_left₀ (by positivity) hInv 2).trans_lt
    (ratio_of_enclosures hD₀.le hD hN)

end TwoShift

#print axioms TwoShift.tail_bound
#print axioms TwoShift.witness_tail_arithmetic
#print axioms TwoShift.constant_of_enclosures
