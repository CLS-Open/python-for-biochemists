# Module 5 — Weighting, error models and the shape of noise

**Time:** 90 minutes. **Prerequisite:** [Module 4](04-uncertainty-standard-errors-and-the-bootstrap.md).

## Learning objectives

By the end of this module you can:

1. state the assumption ordinary least squares makes about the error, and test whether your data
   satisfy it;
2. fit weighted by `1/sd²`, and interpret the reported `scale factor`;
3. compute an LOD and LOQ with an explicitly named σ, and explain why the same standards give
   two different detection limits;
4. explain why a Hill coefficient is a shape parameter rather than evidence of cooperativity.

---

## 5.1 The assumption you never signed up for

Ordinary least squares minimises `Σ(yᵢ − ŷᵢ)²`. Doing that is only statistically efficient
if every point has the **same** error variance. That is an assumption about your assay, and
in enzyme kinetics it is usually false.

Two realistic error shapes:

| Shape | Where it comes from | Consequence for OLS |
|---|---|---|
| Constant (homoscedastic) | detector noise floor dominates | none — OLS is correct |
| Proportional to the signal (heteroscedastic) | pipetting error, dilution error, path-length variation | OLS is dominated by the largest signals |

`datasets/sample/calibration_standards.csv` was generated with 4 % proportional error plus a
small absolute floor, which is the second row of that table. Its replicate SD rises from
about 0.001 near the blank to about 0.4 absorbance units at 200 µM.

## 5.2 What weighting does

Weighted least squares minimises `Σ wᵢ(yᵢ − ŷᵢ)²` with `wᵢ = 1/σᵢ²`. Points whose error you
know to be large stop shouting.

The package exposes it as a `sigma` argument (SciPy's convention: pass standard deviations,
not variances):

```python
import csv
from kinetics.fitting import fit_michaelis_menten

rows = list(csv.DictReader(open("datasets/sample/michaelis_menten_mean.csv", encoding="utf-8")))
S = [float(r["S"]) for r in rows]
v = [float(r["v"]) for r in rows]
sd = [float(r["sd"]) for r in rows]
unweighted = fit_michaelis_menten(S, v)
weighted = fit_michaelis_menten(S, v, sigma=sd)
print(f"unweighted: Vmax {unweighted.params['vmax']:.3f}  Km {unweighted.params['km']:.3f}")
print(f"weighted:   Vmax {weighted.params['vmax']:.3f}  Km {weighted.params['km']:.3f}")
```

```
unweighted: Vmax 9.984  Km 24.877
weighted:   Vmax 9.761  Km 23.922
```

The replicate SDs in that file run from 0.020 to 0.313, so the weights are far from equal
and the weighted fit moves: Vmax by 2 %, Km by 4 %. Against the simulated truth (10.0,
25.0) the *unweighted* fit is closer here — 9.984 against 9.761 — which is the honest
observation and not the moral of the module. On one dataset, every additional decision has
a cost; weighting wins on efficiency across many datasets, given that the weights are
right, and it can lose on any single one.

Where the decision stops being marginal is a calibration over two orders of magnitude.
`datasets/sample/calibration_mean.csv` reports a detection limit 58× larger when
unweighted: LOD 6.36 µM against 0.109 µM (section 5.4 and the
[Figure 4 caption](../index.md)). That is not a 2 % difference; it is the difference
between "this assay can see 0.1 µM" and "this assay can see 6 µM".

## 5.3 Choosing the weights, honestly

In descending order of defensibility:

1. **Measured replicate SDs.** You have them, they are in the file, use them. The package
   takes `sigma=` or, for calibration, `replicate_sd=`.
2. **A documented functional form.** "Weights are 1/v² because the assay's error is
   proportional to the signal" is a statement a reviewer can check. The calibration module
   implements `weighted=True` as `1/signal²` and records that choice in a note.
3. **Inverse variance from a pilot study.** Fine, if you say so.
4. **Equal weights.** The default, and a legitimate choice if the noise really is constant —
   but say *that*, not "we used the standard method".

What is not defensible: choosing weights after looking at which ones make the fit look
better. That is fitting the uncertainty to the residuals.

## 5.4 The calibration case, in full

```python
import numpy as np
from kinetics.calibration import fit_calibration
import csv

rows = list(csv.DictReader(open("datasets/sample/calibration_mean.csv", encoding="utf-8")))
x = np.array([float(r["concentration"]) for r in rows])
y = np.array([float(r["signal"]) for r in rows])
sd = np.array([float(r["sd"]) for r in rows])

for label, kwargs in (("unweighted", {}), ("weighted", {"weighted": True, "replicate_sd": sd})):
    r = fit_calibration(x, y, **kwargs)
    print(f"{label:<11} slope {r.slope:.5f}  LOD {r.lod:8.4f} µM  LOQ {r.loq:8.4f} µM")
```

```
unweighted  slope 0.04764  LOD   6.3578 µM  LOQ  19.2660 µM
weighted    slope 0.04901  LOD   0.1093 µM  LOQ   0.3312 µM
```

The simulated truth is slope 0.0485. The weighted slope is off by 1.1 %; the unweighted by
1.8 %. The interesting number is the detection limit: a factor of 58.

**Why the ratio is so large.** `LOD = 3.3σ/slope`, and the unweighted fit estimates σ from
all residuals pooled, which is dominated by the 200 µM standards where the absolute
scatter is largest (σ ≈ 0.4). That σ is then applied to the *low* end, where the true
noise is 0.001. The weighted fit, which knows the noise at each standard, reports σ at the
lowest standard instead.

Neither number is a "bug": each is the honest consequence of its assumption. The
unweighted one is the answer to "what is my detection limit if I insist the noise is the
same everywhere", and it is 58× wrong about the chemistry.

## 5.5 The trap this module is really about

The weighted branch of `fit_calibration` had a units bug during development: with weights
`1/sd²`, `sqrt(weighted RSS / dof)` is **dimensionless** (it is the ratio of observed to
assumed scatter, ≈ 1.0 when the weights are right). Dividing by the slope produced an LOD
of 66 µM — a plausible-looking number, in the right units, wrong by 600×.

It was caught by a test that asserted the reported σ must be on the scale of the replicate
SDs:

```python
assert 0.001 < res.residual_sd < 0.01
```

The lesson generalises: when a number looks reasonable, check its *units* against something
you can measure. `tests/test_calibration.py::test_weighted_residual_sd_is_in_signal_units_not_a_scale_factor`
is that check, and it is named for the bug.

## 5.6 Hill coefficients: a shape parameter, not a mechanism

The Hill equation `v = Vmax·Sⁿ/(K½ⁿ + Sⁿ)` fits sigmoid data well. Two cautions:

* **n is a shape descriptor.** `n = 2.04 ± 0.11` on `datasets/sample/cooperative_hill.csv`
  means the curve is steeper than a hyperbola. It does **not** demonstrate subunit
  cooperativity: one site with a bad assay, substrate activation, or a mixture of two
  enzymes can produce the same curve.
* **n < 1 is a warning about the assay**, not a discovery of negative cooperativity. It
  usually means substrate inhibition or a rate measured after the initial-rate window.

The package's Hill fit carries this note automatically:

```
a Hill coefficient is a curve-shape parameter until a mechanism is demonstrated
```

And `tests/test_fitting.py` measures what AICc does with the extra parameter: on genuine
Michaelis–Menten data the Hill model wins AICc in 2 of 20 seeds; on genuinely sigmoid data
it wins 20 of 20. Both numbers are asserted, so the claim is reproducible rather than
rhetorical.

## 5.7 Weighting and model comparison interact

AICc assumes Gaussian residuals with constant variance. If you compare models fitted with
*different* weightings, the information criteria are on different scales and the comparison
is invalid — the same problem that makes linearised fits incomparable (Module 3).

Rule: **compare models fitted to the same response with the same weights.** `compare_inhibition_models`
does exactly that: one `_prepare` call, one weighting, three mechanisms.

## 5.8 What to write

> Standards (0–200 µM, n = 3 replicates each) gave a linear response
> (`signal = 0.04901·C + 0.01234`, R² = 0.998). Because the replicate SD increased with
> concentration (0.0016 at the blank to 0.4 A₃₄₀ at 200 µM), the line was fitted by
> weighted least squares with `w = 1/SD²`. LOD and LOQ were computed as `3.3σ/S` and
> `10σ/S` with σ taken at the lowest standard (0.00164 A₃₄₀), giving 0.109 µM and
> 0.331 µM. An unweighted fit of the same data would report an LOD of 6.36 µM; that
> difference is a property of the weighting assumption, not of the assay.

## 5.9 Exercises

1. Rerun the calibration with `weighted=True` but no `replicate_sd`. The module uses
   `1/signal²`; how does the LOD change, and which assumption is now in force?
2. Fit `cooperative_hill.csv` with the Hill equation fixed at `n = 1` and compare AICc
   with the free fit. Which is preferred, and by how much?
3. Add a single wild point (double the value at 10 µM) to the Michaelis–Menten mean data.
   Fit unweighted, then weighted with σ inflated only at that point. Which recovers the
   known truth more closely?
4. Why must you not compare AICc between an unweighted and a weighted fit of the same data?

## 5.10 What this module did **not** cover

* Iteratively reweighted least squares, where the weights are estimated from the fit and
  then updated.
* Robust regression (Huber, Tukey biweight), which limits the influence of outliers rather
  than down-weighting them by a known σ.
* Measurement-error models where the x variable also has error (errors-in-variables); the
  package assumes x is exact, which is right for `[S]` and wrong for anything you measured
  with error.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| The weighted fit is worse than the unweighted one | The declared `sd` values do not describe the actual scatter (check the scale factor) | Recompute the SDs from the raw replicates, or fit unweighted |
| `LOD` larger than the top standard | σ is a residual SD dominated by high-concentration points | Use the SD at the lowest standard, and say which σ you used |
| `LOQ` smaller than `LOD` | You mixed multipliers (3.3 and 10) between the two calculations | Both come from the same σ; check the arithmetic |
| Every weight is identical | The `sd` column is constant, so weighting changes nothing | That is not weighting; say the fit is unweighted |

---

**Next:** [Module 6 — choosing between models without fooling yourself](06-choosing-between-models-without-fooling-yourself.md)
