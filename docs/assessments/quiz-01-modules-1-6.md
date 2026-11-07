# Quiz 1 — Modules 1–6

**Time:** 60 minutes. **Closed book** except for one page of notes. **Total: 50 marks.**

Answer in sentences, not keywords. Every numerical answer should carry its unit, and every
"which is better" answer must say *better at what*.

Run nothing for questions 1–6; questions 7–12 ask you to predict output you have seen, which
is a different skill from reproducing it.

---

## Section A — What the numbers mean (12 marks)

**1.** (3 marks) A dataset records `v = 7.97 µmol min⁻¹ (mg protein)⁻¹` at `[S] = 100 µM`.
State (a) the reaction rate in mol s⁻¹ per mg of protein, (b) one thing you cannot tell from
this number alone, and (c) why the denominator matters when comparing two published Km values.

**2.** (3 marks) Vmax is described as *an extrapolation, not a measurement*. Give the
`michaelis_menten.csv` fit's Vmax and Km from the worked example, then state what substrate
concentration would be needed to measure `v = 0.99·Vmax` directly, in terms of Km.

**3.** (2 marks) Why does `sd` in `michaelis_menten_mean.csv` fail to increase monotonically
with `[S]` (100 µM has 0.31285, 200 µM has 0.25232)?

**4.** (4 marks) The `kinetics mm` output for `michaelis_menten_mean.csv` prints
`weighting: 1/sd² from the sd column`. Explain in one sentence what that choice assumes, and
give one reason a reader could not reproduce the fit if the report had omitted that row.

## Section B — Transformations (10 marks)

**5.** (4 marks) On `michaelis_menten_mean.csv`, `kinetics compare` returns Vmax = 9.98381 for
the nonlinear fit and 11.9869 for Lineweaver–Burk, against a simulated truth of 10.0. Compute
the relative error of each. Which of the three linearisations was closest, and by how much?

**6.** (3 marks) At `[S] = 1 µM`, `v = 0.3601` with `sd = 0.02161`. Show that propagating the
error through `y = 1/v` gives σ ≈ 0.1666, and state what fraction of `1/v` that is.

**7.** (3 marks) A colleague says: "I ran the double reciprocal on my data and the line was
beautiful, R² = 0.999, so the bias isn't a problem here." Using the repository's own simulation
result (how often Lineweaver–Burk recovers Vmax within 1 % across 200 simulated assays), give
the two-sentence reply.

## Section C — Uncertainty (10 marks)

**8.** (3 marks) `fit_michaelis_menten` on the 27-row series gives Vmax = 9.9838 ± 0.1069. The
percentile bootstrap (2000 resamples, seed 20261107) gives 9.7095 – 10.2630. The normal-theory
interval is 9.774 – 10.193. State which you would report and why, referring to the shape of
the interval.

**9.** (3 marks) Fitting the 9 means gives an identical point estimate to four decimal places
but a *wider* bootstrap interval (9.6349 – 10.3246). Explain why more rows produced a narrower
interval while the same 9 means produced a wider one.

**10.** (2 marks) Give one source of error that a bootstrap of these rows cannot represent.

**11.** (2 marks) Why does `bootstrap_ci` make `seed` a required argument? Give the sentence you
would write in a methods section that satisfies it.

## Section D — Choosing and refusing (18 marks)

**12.** (4 marks) Two models are fitted to the same response: competitive inhibition with
AICc = −83.17 and noncompetitive with AICc = −40.62. State the ΔAICc and the conclusion, in the
form of a sentence that could go in a paper.

**13.** (4 marks) On `inhibition_weak.csv` the ranking is uncompetitive (−54.09),
noncompetitive (−53.62), competitive (−52.41), while the true mechanism is competitive. What do
you report, and what must you *not* report?

**14.** (4 marks) The same seven calibration standards give LOD = 0.109296 weighted and
LOD = 6.357793 unweighted. Explain where the factor of 58 comes from, naming the σ used in each
case.

**15.** (3 marks) A four-parameter mixed-inhibition fit on `inhibition_series.csv` returns
`alpha = 97030109.768` with a standard error of `nan` and AICc −80.01, worse than the
three-parameter competitive fit at −83.17. What does each of those three facts tell you?

**16.** (3 marks) In one sentence each: what is the difference between an LOD and an LOQ, and
what is the difference between "the best model by AICc" and "the mechanism"?

---

# Answer key

**1.** (a) 7.97 µmol min⁻¹ mg⁻¹ = 7.97 × 10⁻⁶ mol / 60 s per mg = **1.33 × 10⁻⁷ mol s⁻¹ mg⁻¹**.
(b) Anything about the enzyme's kinetics — one point on a curve says nothing about Vmax, Km, or
whether the model is Michaelis–Menten. (c) Rates per mg total protein and per mg *pure enzyme*
differ by the purification factor, so two Km values are comparable only if the denominators are;
Km itself is a concentration and does not depend on it, but Vmax does.

**2.** Vmax = **9.9838**, Km = **24.8774 µM** (27-row series). `v = 0.99·Vmax` needs
`S/(Km+S) = 0.99`, i.e. **S = 99·Km ≈ 2,460 µM** — 99 times Km, which is why nobody measures
Vmax directly.

**3.** Two reasons: the generator combines a 3 % proportional term with an absolute floor
(`error_model: proportional_plus_absolute`), and each `sd` is a sample SD of only three
replicates, so it is a noisy estimate. The trend is upward; the sequence is not monotone.

**4.** It assumes the replicate SD is the true measurement uncertainty at each concentration and
that the errors are independent, so each point should contribute in proportion to `1/σ²`. A
reader without that row cannot tell whether the fit minimised weighted or unweighted residuals
— and on this dataset those give different answers (9.76057/23.9224 against 9.98381/24.8774).

**5.** Nonlinear: (9.98381 − 10)/10 = **−0.16 %**. Lineweaver–Burk: (11.9869 − 10)/10 =
**+19.87 %**. **Hanes–Woolf** was closest: 9.98582 → **−0.14 %** (Km 25.0091, +0.04 %).

**6.** `1/v = 2.77701`. σ(1/v) ≈ σ_v/v² = 0.02161/0.3601² = **0.16665**, which is **6.0 %** of
2.77701. At 200 µM the same calculation gives 0.00321 on 0.11276 = 2.8 % — the low point's
transformed error is more than twice as large.

**7.** Across 200 simulated assays with the repository's error model, Lineweaver–Burk recovered
Vmax within 1 % in only **5 cases (2.5 %)**. A beautiful line is what the transformation does
whether or not the parameters are right, so "it looked fine" is a statement about the plot, not
about the estimate.

**8.** The **bootstrap** interval, 9.7095 – 10.2630. It is wider than the normal-theory
interval on both sides (0.2743 below the estimate against a half-width of 0.2095), because Vmax
is a ratio whose sampling distribution is not symmetric or thin-tailed in the way `estimate ±
1.96·SE` assumes. When the two disagree, the bootstrap is the one you can defend.

**9.** The 27 rows contain only 9 conditions; the replicates estimate within-condition noise,
which the 9 means already capture, so the extra rows add little information about where the
curve bends. Resampling 27 rows treats them as 27 independent observations (pseudoreplication)
and therefore *understates* the interval — 0.554 against 0.690 for the means. Tripling the rows
narrowed the interval by 20 %, not the 42 % that independence predicts.

**10.** Any error shared by all replicates: a miscalibrated pipette, a drifting temperature, a
mis-stated stock concentration, or model misspecification. It resamples with itself, so it never
appears in the interval.

**11.** Because a Monte-Carlo interval that cannot be reproduced is not a result — five seeds
move the Vmax endpoints by ~0.03 units. Methods sentence: *"95 % confidence intervals are
percentile bootstrap intervals from 2000 resamples of the substrate–rate pairs, seed 20261107."*

**12.** ΔAICc = **42.55** (competitive −83.17 vs noncompetitive −40.62). Sentence: *"The
competitive model is better supported than the noncompetitive model by ΔAICc = 42.6, which is
far outside the range where the data could be said to be indifferent between them."*

**13.** Report: *"all three mechanisms fit the data equally well (ΔAICc ≤ 1.68); the experiment
cannot distinguish them."* Do **not** report "uncompetitive inhibition" as the answer, and do not
report the Ki values (47.4696 ± 34.1 etc.), which are unidentified.

**14.** Both use `LOD = 3.3σ/S`, so the ratio is the ratio of the σ values: the weighted fit uses
σ = **0.001641**, the SD at the lowest standard (the blank, in signal units), while the
unweighted fit uses the residual SD of the regression, **0.0917743**, which is dominated by the
high-concentration points. 0.0917743/0.001641 = 55.9, and the small difference from 58.17 comes
from the slope (0.0490126 against 0.0476353).

**15.** `alpha → 9.7 × 10⁷` means the optimiser has run to the **competitive limit** (α → ∞), so
the four-parameter model has rediscovered the three-parameter answer. `nan` standard errors mean
the covariance matrix was not estimable — an estimate with no interval is not a measurement.
AICc −80.01 against −83.17 means the extra parameter **cost more than it bought**.

**16.** LOD is where a signal can just be distinguished from the blank (3.3σ); LOQ is where the
relative uncertainty is fit for a measurement (10σ) — on `calibration_mean.csv`, 0.109296
against 0.331201. "Best model by AICc" is a statement about relative fit *among the models
tried*; "the mechanism" is a claim about the enzyme, which no fit can establish on its own.

---

## Marking scheme

| Section | Marks | What earns full credit |
|---|---:|---|
| A — Units and parameters | 12 | Units carried through; Vmax called an extrapolation; the `sd` non-monotonicity explained by the error model *and* the n = 3 noise |
| B — Transformations | 10 | Relative errors computed, Hanes identified, error propagation shown, the simulation number quoted |
| C — Uncertainty | 10 | Asymmetry explained; pseudoreplication named; one genuinely unrepresentable error given |
| D — Choosing and refusing | 18 | ΔAICc interpreted with the right threshold; the non-claim stated for `inhibition_weak`; σ definitions named for the LOD; α read as the competitive limit |

Deduct 1 mark per answer that quotes a number without its unit where a unit exists. Answers
that name a mechanism for `inhibition_weak.csv` score 0 for question 13 regardless of the rest
of the answer.
