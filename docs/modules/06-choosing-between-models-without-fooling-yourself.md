# Module 6 — Choosing between models without fooling yourself

**Time:** 105 minutes. **Prerequisite:** [Module 5](05-weighting-error-models-and-the-shape-of-noise.md).

## Learning objectives

By the end of this module you can:

1. compute and interpret ΔAICc, including the small-sample correction, and say when two models
   are indistinguishable;
2. explain why a linearised fit's information criterion is not comparable with a nonlinear one;
3. type an inhibition mechanism when the design permits it, and refuse when it does not;
4. write the non-claim — the sentence that reports a number and the doubt together.

---

## 6.1 The problem in one sentence

Adding a parameter always improves the fit. The question is whether it improved it *enough*
to justify itself, and no amount of R² will answer that.

## 6.2 AICc, and only as much theory as you need

```
AIC  = n·ln(RSS/n) + 2k
AICc = AIC + 2k(k+1)/(n − k − 1)
```

Read it as **fit minus penalty**: the first term gets smaller as the model explains more;
the second gets bigger as you add parameters. Lower AICc = better trade-off. The correction
term matters whenever `n/k < 40`, which in enzymology is nearly always — with nine data
points and three parameters the penalty is `2·3·4/5 = 4.8` units, larger than the 2 units
AIC would charge.

Two things AICc is **not**:

* a probability that the model is true (it is a relative measure among the candidates you
  supplied);
* comparable across different responses, different weightings, or transformed data.

The package raises rather than returning a number when the correction is undefined
(`n − k − 1 ≤ 0`):

```python
from kinetics.fitting import aicc

try:
    aicc(4.0, 3, 2)
except ValueError as exc:
    print("refused:", exc)
```

```
refused: AICc undefined: n=3, k=2 leaves no degrees of freedom
```

## 6.3 The Δ scale, with the numbers this repository measured

Conventional reading (Burnham & Anderson 2004): Δ ≤ 2 → the data do not separate the
models; 4–7 → considerably less support; > 10 → effectively none.

Two measurements from this repository, both asserted in tests:

| Comparison | ΔAICc | Verdict |
|---|---:|---|
| Hyperbola vs Hill on `cooperative_hill.csv` (true n = 2) | 26.97 | Hill, decisively |
| Competitive vs non-competitive on `inhibition_series.csv` (true competitive) | 42.55 | competitive, decisively |
| Competitive vs non-competitive on `inhibition_weak.csv` (true competitive) | 1.68 | nothing is decided |

The third row is the important one. Same mechanism, same generator, different experimental
design, and the honest answer changes from "competitive" to "unidentified".

## 6.4 Fitting the mechanism you cannot see: inhibition typing

The (S, I, v) surface is fitted with each candidate mechanism and the fits are compared:

* **competitive** — inhibitor binds free enzyme; apparent Km rises, Vmax unchanged;
* **non-competitive** — binds E and ES equally; Vmax falls, Km unchanged;
* **uncompetitive** — binds ES only; both apparent Km and Vmax fall.

```bash
kinetics inhibition datasets/sample/inhibition_series.csv
```

```
| best by AICc | competitive |
| caveat | a best-fitting mechanism is not a proven mechanism |

| model | method | n | k | R² | AICc |
|---|---|---:|---:|---:|---:|
| competitive | nonlinear | 20 | 3 | 0.99811 | -83.17 |
| noncompetitive | nonlinear | 20 | 3 | 0.98412 | -40.62 |
```

The fitted parameters are Vmax 9.61, Km 24.04, Ki 8.59 against simulated truth 10.0, 25.0,
8.0 — all within 4 %, which is what a well-designed experiment should deliver.

Now the same command on the badly designed experiment:

```bash
kinetics inhibition datasets/sample/inhibition_weak.csv
```

```
| best by AICc | uncompetitive |
```

The best model is **wrong**. It is uncompetitive, while the data were generated with a
competitive inhibitor (recorded in `datasets/metadata/inhibition_weak.yml`). ΔAICc for
competitive is 1.68, for non-competitive 0.48 — all three inside the "indistinguishable"
band, and the notes say so:

> indistinguishable models (ΔAICc ≤ 2): competitive, noncompetitive — report the ambiguity
> rather than the winner

Why the design failed: the inhibitor concentrations used (0–2 µM) never exceed 0.004·Ki
(Ki = 500 µM), so the mechanism barely perturbs the rate and the three models are nearly
the same surface. **No amount of fitting recovers information the experiment did not
collect.** That is not a software limitation; it is the definition of underpowered.

## 6.5 The general form of the trap

The weak-inhibition dataset is a template for a class of mistakes:

| Situation | What a careless report says | What the data support |
|---|---|---|
| Two models within ΔAICc 2 | "Model A fits best" | "the data do not separate A and B" |
| Hill n = 1.15 ± 0.40 | "positive cooperativity" | "no evidence of deviation from a hyperbola" |
| Mixed inhibition with α = 0.3 ± 4 | "mixed inhibition with α = 0.3" | "α is unidentified on this data" |
| R² = 0.997 | "excellent linearity" | "the residuals are small compared to the range, nothing more" |

The package helps with all four by putting a note in the result next to the numbers, but it
cannot stop you writing the careless version in a report. That is the human half of the
skill.

## 6.6 Over-parameterisation, concretely

`mixed_inhibition` has four parameters (Vmax, Km, Ki, α). Fitting α needs enough inhibitor
concentrations on both sides of Ki to see whether the Vmax asymptote moves. With the
`inhibition_series` design (four [I] from 0 to 2·Ki) you can type the mechanism; you cannot
estimate α with any precision. The module does not offer a mixed-inhibition fit for exactly
this reason — offering it would invite a number with a 95 % interval spanning two orders of
magnitude to be reported as a result.

If you need α, the package gives you the forward model (`mixed_inhibition`) and
`fit_nonlinear`, and the responsibility.

```python
import csv
from kinetics.fitting import fit_nonlinear
from kinetics.models import mixed_inhibition

rows = list(csv.DictReader(open("datasets/sample/inhibition_series.csv", encoding="utf-8")))
S = [float(r["S"]) for r in rows]
I = [float(r["I"]) for r in rows]
v = [float(r["v"]) for r in rows]


def surface(s, vmax, km, ki, alpha, _I=I):
    return mixed_inhibition(s, vmax, km, _I, ki, alpha)


fit = fit_nonlinear(
    S, v, surface, (10.0, 25.0, 8.0, 2.0), ("vmax", "km", "ki", "alpha"), name="mixed_inhibition"
)
for name, value in fit.params.items():
    print(f"{name:<6} {value:8.3f} ± {fit.stderr[name]:7.3f}")
print("AICc", round(fit.aicc, 2))
```

The output of that block, run on the committed data:

```
vmax     9.606 ±     nan
km      24.037 ±     nan
ki       8.591 ±     nan
alpha  97030109.768 ±     nan
AICc -80.01
```

Read it carefully, because it is the clearest negative result in the repository:

* `alpha` has run away to 9.7 × 10⁷. Since α → ∞ is exactly the competitive limit, the
  optimiser has told you that the data prefer "no binding to ES" — it has rediscovered the
  competitive mechanism using four parameters instead of three.
* Every standard error is `nan`, and the result carries the note *"the optimiser could not
  estimate the covariance matrix, so every standard error in this result is NaN"*. An
  estimate with no interval is not a measurement.
* AICc is −80.01, worse than the three-parameter competitive fit's −83.17: the extra
  parameter costs more than it buys, which is the argument of this whole module.

The honest sentence is not "α = 0.3, suggesting mixed inhibition". It is "the data are
consistent with the competitive limit and do not constrain α", followed by a recommendation
to use the three-parameter model.

## 6.7 What to report when you cannot decide

A perfectly respectable results paragraph:

> The (S, I, v) surface (5 substrate × 3 inhibitor concentrations, 15 points) was fitted
> with competitive, non-competitive and uncompetitive mechanisms. All three fitted
> adequately (R² > 0.997) and the smallest ΔAICc was 0.48 (uncompetitive vs
> non-competitive), well inside the range where the data do not discriminate. The
> experiment therefore does not identify the mechanism; a design with inhibitor
> concentrations spanning Ki would be required.

This is a *better* paragraph than one that names a winner, because it tells the next person
what to do.

## 6.8 A checklist before you choose a model

1. Are the candidate models fitted to the same response, the same weights, the same points?
2. Is AICc defined for every candidate (n − k − 1 > 0)?
3. Is any ΔAICc ≤ 2? If so, report the ambiguity.
4. Does the winning model's extra parameter have an interpretable estimate with an interval
   narrower than "anything"?
5. Would a different plausible candidate give the same prediction over your data range? If
   yes, the choice is not supported by these data.

## 6.9 Exercises

1. Compute ΔAICc between the Hill fit and the Michaelis–Menten fit on
   `cooperative_hill.csv` from the printed values (AICc −22.59 and +4.37). Which columns of
   the CLI output would you quote?
2. Simulate `inhibition_series` with [I] = 0, 0.5, 1, 2 µM instead of 0, 2, 8, 16 µM. Run
   `compare_inhibition_models`. Where does the comparison stop working?
3. Take the `inhibition_weak` result and write the two-sentence "cannot decide" version of
   the conclusion. Then write the careless version. Note what information the careless one
   destroys.
4. Why is a linearised fit excluded from the CLI's AICc ranking? Give the units argument.

## 6.10 What this module did **not** cover

* Cross-validation and information criteria other than AICc (BIC, WAIC, LOO).
* Model averaging, which is the principled response to ΔAICc ≤ 2 in prediction tasks.
* F-tests for nested models, which answer a narrower question and are widely misused.
* Mechanism discrimination by isotope effects, dead-end inhibitors, or product-inhibition
  patterns — the experiments that actually answer the question this module can only flag as
  unidentified.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| AICc reads "not comparable" | One fit is on a transformed response (1/v, v/S, S/v) | Compare like with like; report linearised R² separately |
| `AICc undefined: n=3, k=2 leaves no degrees of freedom` | Too few points for the correction | Report no comparison; that is the honest answer |
| Two mechanisms tie within Δ = 0.5 | The design cannot separate them, usually too few [I] levels | Report the tie; design a better experiment (Practical 5) |
| Every standard error is `nan` in a four-parameter fit | The covariance matrix was not estimable — the model is over-parameterised | Drop to three parameters; read Module 6 §6.6 |

---

**Next:** [Practical 5 — typing an inhibitor](../practicals/practical-05-typing-an-inhibitor.md) or the
[quiz](../assessments/quiz-01-modules-1-6.md).
