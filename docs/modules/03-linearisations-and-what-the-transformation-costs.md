# Module 3 — Linearisations and what the transformation costs

**Time:** 90 minutes. **Prerequisite:** [Module 2](02-the-michaelis-menten-model-and-the-two-numbers-people-misquote.md).

## Learning objectives

By the end of this module you can:

1. write the Lineweaver–Burk, Hanes–Woolf and Eadie–Hofstee transformations and say which axis
   carries the error in each;
2. propagate the error through `1/v` and state the fraction it represents at two different
   substrate concentrations;
3. measure the bias instead of asserting it, using the 400-assay distribution and the 200-assay
   failure count;
4. name a legitimate use of a linearisation — as a diagnostic — and the conditions under which
   a Hanes–Woolf fit is an acceptable approximation.

---

## 3.1 Three ways to turn a curve into a line

Before computers, a straight line drawn by eye was the only practical way to get Vmax and
Km, so the hyperbola was rearranged into linear forms:

| Name | Axes | Slope | Intercept | Reference |
|---|---|---|---|---|
| Lineweaver–Burk | `1/v` vs `1/S` | `Km/Vmax` | `1/Vmax` | Lineweaver & Burk 1934 |
| Hanes–Woolf | `S/v` vs `S` | `1/Vmax` | `Km/Vmax` | Hanes 1932 |
| Eadie–Hofstee | `v` vs `v/S` | `−Km` | `Vmax` | Eadie 1942; Hofstee 1959 |

All three are algebraically exact rearrangements of `v = Vmax·S/(Km + S)`. If the data had
no error, all three would give the same answer, exactly. The differences are entirely
about **what happens to the error**.

## 3.2 The error propagation, in one line each

Ordinary least squares minimises the sum of squared *vertical* distances. Apply a
transformation `y' = f(y)` and you are now minimising distances in `y'`, which is not the
same as minimising them in `y`.

* **Lineweaver–Burk**: `y' = 1/v`, so `dy'/dy = −1/v²`. An error of `σ` in `v` becomes
  `σ/v²` in `1/v`. At `S = 1 µM` and `v = 0.36`, `1/v = 2.78` and an absolute error of 0.05
  becomes 0.39 in the transformed value — an eight-fold amplification, applied to the points
  that were *already* the least reliable (small `v`, small `S`).
  Then the regression reads `1/Vmax` off an intercept obtained by extrapolating back to
  `1/S = 0`, i.e. to `S = ∞`, from points whose transformed errors are largest exactly where
  the extrapolation is anchored.
* **Hanes–Woolf**: `y' = S/v`, so `dy'/dv = −S/v²`. The error still grows as `v` falls, but
  the x-axis is `S`, which is the variable you set and therefore the one without
  measurement error. This is the least damaging of the three.
* **Eadie–Hofstee**: no reciprocal of the response, but `v` appears on *both* axes
  (`v` against `v/S`), so the errors in the x and y variables are correlated, which violates
  an assumption of ordinary least squares. It bends the estimate towards low `v`.

## 3.3 The claim, measured

Documentation that says "Lineweaver–Burk is biased" is not evidence. `tests/test_linearisation.py`
runs 300 seeded assays with **constant-variance** noise at σ = 0.1 (1 % of Vmax) and
records every estimate. The measured result:

| Method | median \|ΔVmax\| | median \|ΔKm\| | within 20 % of Vmax | catastrophic (>100 % error) |
|---|---:|---:|---:|---:|
| nonlinear | 0.121 | 0.840 | 100 % | 0 % |
| Hanes–Woolf | 0.206 | 1.695 | 100 % | 0 % |
| Eadie–Hofstee | 0.530 | 2.812 | 97 % | 0 % |
| Lineweaver–Burk | 3.980 | 13.130 | 22 % | 30 % |

Those numbers are asserted in the test file, so if a change to the code makes them untrue
the suite fails. `scripts/make_figures.py` produces the same experiment at 400 seeds with a
different error model (3 % proportional plus 0.05 absolute) for Figure 2: there the median
Lineweaver–Burk estimate is 9.74 against a truth of 10.0 — close — but only 57 % of
estimates land within 20 %, and the distribution has a tail that reaches past 25.

Both are the same finding seen from different angles: **the bias is not a shift in the
median on every dataset; it is the variance and the tails.**

## 3.4 What it looks like on a single dataset

`datasets/sample/michaelis_menten_mean.csv` is one well-behaved dataset (nine
concentrations, Km ≈ 25 µM, highest `[S]` = 8·Km):

```python
import csv
from kinetics.fitting import compare_linearisations

rows = list(csv.DictReader(open("datasets/sample/michaelis_menten_mean.csv", encoding="utf-8")))
S = [float(r["S"]) for r in rows]
v = [float(r["v"]) for r in rows]
for result in compare_linearisations(S, v):
    print(f"{result.model:<18} Vmax {result.params['vmax']:6.3f}   Km {result.params['km']:6.3f}")
```

```
michaelis_menten   Vmax  9.984   Km 24.877
hanes              Vmax  9.986   Km 25.009
eadie_hofstee      Vmax 10.013   Km 25.163
lineweaver_burk    Vmax 11.987   Km 31.905
```

Every transform agrees to within 0.2 on Vmax except Lineweaver–Burk, which is 20 % high —
on a dataset where the nonlinear fit is 0.16 % from the simulated truth. The order of
agreement is exactly the order the error propagation predicts.

## 3.5 Why the R² looks great anyway

Run the CLI on the same file:

```bash
kinetics compare datasets/sample/michaelis_menten_mean.csv
```

The comparison table shows R² values of 0.99934 (nonlinear), 0.99970 (Hanes),
0.99140 (Eadie–Hofstee) and 0.99875 (Lineweaver–Burk). The *worst* method for Vmax has a
higher R² than the best one, because R² is computed on the transformed scale where the
regression was performed. A high R² in a Lineweaver–Burk paper means "the reciprocals made
a straight line", which is guaranteed by the algebra, not by the quality of the estimate.

You will also see the AICc column says **not comparable** for the linearised fits. That is
deliberate: RSS on `1/v` and RSS on `v` are in different units, so information criteria
cannot be compared across them. The CLI excludes these fits from its AICc ranking and says
so; `docs/validation.md` explains what happens if you ignore that.

## 3.6 When is a linearisation acceptable?

Never as the final estimate if you have a computer. Acceptable as:

* **A diagnostic.** Eadie–Hofstee is still the best plot for *seeing* deviations
  (two lines for two enzymes, curvature for cooperativity, a bend at high `S` for substrate
  inhibition) because both axes are in rate units and the eye is good at lines.
* **A historical comparison.** If you must compare your Km with a 1970s value, that value
  came from a double-reciprocal plot, so its uncertainty is not the ± you see printed.
* **A teaching device.** This module exists because the transformation is a compressed
  lesson in error propagation.

The package keeps all three fits — `fit_lineweaver_burk`, `fit_eadie_hofstee`,
`fit_hanes` — and each result carries a note describing its bias, because deleting the
methods would delete the lesson.

## 3.7 The degenerate case that will bite you

Lineweaver–Burk on exactly proportional data (`v/S` constant) has an intercept of zero, so
`Vmax = 1/intercept = ∞`:

```python
from kinetics.fitting import fit_lineweaver_burk

result = fit_lineweaver_burk([1, 2, 5, 10, 20], [2, 4, 10, 20, 40])
print(result.params["vmax"], result.params["km"])
print(result.notes[-1][:120])
```

```
inf nan
the fitted intercept is zero to machine precision, so Vmax is unbounded and no Km
```

The keys are kept (so downstream code does not raise `KeyError`) and the values are
`inf`/`nan` rather than a made-up number. This is the shape of every "why is my Vmax
1e17?" question.

## 3.8 What to do instead

1. Fit the hyperbola directly (`fit_michaelis_menten`) — it is one function call.
2. If you must use a linearisation for a plot, fit the nonlinear model for the numbers.
3. Report the linearisation's bias note if you show it. The package's results carry the
   note automatically; a human author should carry it too.

## 3.9 Exercises

1. Reproduce the 300-assay experiment in `tests/test_linearisation.py` with σ = 0.05 and
   σ = 0.2. Does the ratio of Lineweaver–Burk failures rise faster than the noise?
2. Show algebraically that `1/v` has variance `σ²/v⁴` when `v` has variance `σ²`. Which
   points get the largest weight in the regression, and why is that the wrong way round?
3. Construct a dataset where Lineweaver–Burk *underestimates* Vmax. (Hint: put the noise
   where it hurts.) What does that tell you about "correcting" for the bias by multiplying
   by a constant?
4. On the same dataset, Hanes–Woolf and the nonlinear fit agree to 0.002. Does that make
   Hanes–Woolf a safe default?

   *(No: it means this dataset is favourable. Figure 2's 400 assays show Hanes–Woolf with
   an IQR about three times the nonlinear fit's, and Module 5 shows that the choice
   interacts with the error model.)*

## 3.10 What this module did **not** cover

* Weighted linearisations (there is a literature; it is a patch on a bad idea).
* Direct linear plot (Cornish-Bowden & Eisenthal 1974), which is a genuinely different and
  better graphical approach — mentioned in [references](../references.md), not implemented.
* Robust regression, which would reduce the tails but change the estimand.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Lineweaver–Burk returns `vmax = inf, km = nan` | The data are proportional (`v/S` constant), so the double reciprocal has no intercept | That is the correct answer to a degenerate dataset; it is not a bug |
| `ValueError: ... requires positive [S] and positive rates` | The transformation reciprocates, so a zero is fatal | Drop the point *and say you did*, or use the nonlinear fit |
| Hanes–Woolf refuses with "all rates must be positive" | It divides by `v` | Same choice as above |
| `ValueError: ... slope is not identifiable` | Every x value is identical — a blocked experiment | Check the transformation you chose |

---

**Next:** [Module 4 — uncertainty, standard errors and the bootstrap](04-uncertainty-standard-errors-and-the-bootstrap.md)
