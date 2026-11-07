# Module 4 — Uncertainty: standard errors, replicates and the bootstrap

**Time:** 90 minutes. **Prerequisite:** [Module 3](03-linearisations-and-what-the-transformation-costs.md).

## Learning objectives

By the end of this module you can:

1. explain what an asymptotic standard error assumes, and why it is optimistic for a ratio like
   Vmax/Km;
2. run a percentile bootstrap, report its interval, and state the seed as part of the result;
3. distinguish replicate variability from condition variability, and explain why 27 rows can
   carry less information than they appear to;
4. say what a bootstrap cannot fix — shared systematic error and model misspecification.

---

## 4.1 "The Km is 24.9" is not a result

It is an estimate reported without its uncertainty, and it invites a reader to believe
something the data cannot support. Every number in this course that came from a fit comes
with one of three kinds of interval, and they are not interchangeable:

| Kind | What it estimates | Assumes | Cost |
|---|---|---|---|
| Asymptotic standard error | movement of the estimate under repeated sampling, *if the model is right* | correct model, normal-ish errors, large n | free with the fit |
| Bootstrap interval | same, estimated from your own data by resampling | your sample is representative | a second of CPU, a seed |
| Replicate-based interval | spread between independent experiments | you actually did independent experiments | laboratory time |

## 4.2 Asymptotic standard errors, and why they are optimistic

Nonlinear least squares gives standard errors from the curvature of the residual sum of
squares at the minimum: a sharp minimum means a well-determined parameter. Three ways this
misleads:

1. **It assumes the model is right.** Fit a hyperbola to sigmoid data and you get a tight
   standard error around a wrong Km (Module 2 measured a Km inflated by 77 % with
   R² = 0.95).
2. **It is symmetric.** Km is a ratio-like quantity; its sampling distribution is skewed.
   A symmetric ±1.2 around 24.9 implies that a Km of 23.7 is as likely as 26.1, which is
   not true.
3. **It is optimistic at small n.** The formula uses a linear approximation that is worst
   when the data are few — which is the normal situation in enzymology.

The package returns them anyway, because they are the standard and a reviewer will ask for
them:

```python
import csv
from kinetics.fitting import fit_michaelis_menten

rows = list(csv.DictReader(open("datasets/sample/michaelis_menten.csv", encoding="utf-8")))
S = [float(r["S"]) for r in rows]
v = [float(r["v"]) for r in rows]
fit = fit_michaelis_menten(S, v)
print(f"Vmax {fit.params['vmax']:.3f} ± {fit.stderr['vmax']:.3f}")
print(f"Km   {fit.params['km']:.3f} ± {fit.stderr['km']:.3f}")
```

```
Vmax 9.984 ± 0.107
Km   24.877 ± 0.840
```

## 4.3 The bootstrap, in one paragraph and one command

Resample the `(S, v)` pairs with replacement, refit, repeat a few thousand times, and take
the 2.5th and 97.5th percentiles of the resulting estimates. That is the whole method.

```python
from kinetics.fitting import bootstrap_ci

ci = bootstrap_ci(S, v, n_resamples=2000, seed=20261107)
print("Vmax", tuple(round(x, 3) for x in ci["vmax"]))
print("Km  ", tuple(round(x, 3) for x in ci["km"]))
```

```
Vmax (9.71, 10.263)
Km   (23.359, 26.449)
```

Compare with the asymptotic intervals: `9.984 ± 0.107` is roughly `9.77–10.20`, and the
bootstrap gives 9.71–10.26 — slightly wider, and asymmetric around the point estimate
(0.27 below, 0.28 above for Vmax; 1.52 below and 1.57 above for Km). The asymmetry is the
honest part.

## 4.4 The seed is part of the result

`bootstrap_ci` requires a seed and refuses without one. This is not bureaucracy:

* a percentile from 2000 resamples has Monte-Carlo error of its own — rerun with a
  different seed and the interval moves in the second decimal;
* an interval that cannot be reproduced is not a result, it is a memory;
* a reported seed costs one line in a methods section.

`tests/test_fitting.py` asserts that the same seed gives exactly the same interval and that
a different seed gives a different one of similar width. That is the contract.

## 4.5 What the bootstrap does *not* fix

This is the part learners skip and reviewers catch.

| Problem | Does the bootstrap fix it? |
|---|---|
| Skewed parameter distribution | **Yes** — this is what it is for |
| Non-normal residuals | Partly: it resamples your residuals, whatever they are |
| A wrong model | **No.** Bootstrapping a wrong fit gives a confident interval around a wrong number |
| Systematic error (wrong blank, degraded substrate) | **No** |
| Pseudoreplication (replicates that are not independent) | **No** — it will *shrink* the interval, which is worse than useless |
| Small n | **No.** Four points resampled 2000 times is four points |

The pseudoreplication point deserves a worked example. `datasets/sample/michaelis_menten.csv`
contains three replicates at nine concentrations. Treating all 27 rows as independent
observations, as the snippet above did, is *not* the same as having 27 independent
experiments: the three replicates at one concentration share a pipetting error, a
temperature drift, and a cuvette. The bootstrap over 27 rows gives ±0.107 on Vmax; the
bootstrap over the nine *concentrations* would be wider.

The package cannot detect this — it sees rows. You have to know what a row is. (This is the
same error Hurlbert 1984 named; see [references](../references.md).)

## 4.6 Replicates: the interval you can actually defend

If you have three independent repeats of the whole experiment, fit each separately and look
at the spread of the three fitted parameters. That is a small, honest interval that no
statistical machinery can inflate or deflate.

```python
import numpy as np

fits = []
for seed in (1, 2, 3):
    rng = np.random.default_rng(seed)
    noisy = np.array(v) * (1 + 0.03 * rng.standard_normal(len(v)))
    fits.append(fit_michaelis_menten(S, noisy).params["km"])
print("Km over three simulated repeats:", [round(x, 2) for x in fits])
```

Three numbers are enough to see the spread and not enough to estimate a confidence level.
Report the range, say n = 3, and do not attach "95 %" to it.

## 4.7 What to write

> Km was estimated at 24.9 µM (asymptotic SE 0.84, n = 9 concentrations × 3 replicates).
> A percentile bootstrap (2000 resamples, seed 20261107) gave a 95 % interval of
> 23.4–26.4 µM. Both intervals treat the nine concentrations as independent and the
> replicates as measurement error; if the replicates share a systematic error, the true
> uncertainty is larger than either.

That last sentence is the one that distinguishes a careful report from a careless one.

## 4.8 Exercise: make the interval lie

```python
from kinetics.fitting import bootstrap_ci

# constant response: nothing to identify, but the estimator still returns something
try:
    print(bootstrap_ci([1, 2, 3, 4], [5, 5, 5, 5], n_resamples=100, seed=1))
except ValueError as exc:
    print("refused:", exc)
```

```
refused: every bootstrap resample failed to fit; the data do not identify the model
```

The library refuses when the resamples cannot be fitted. A version that silently returned
the point estimate would produce a beautiful zero-width interval for a failed experiment.

## 4.9 Exercises

1. Increase `n_resamples` from 200 to 2000 to 20000 on the same seed. How much does the
   Km interval move? At what point does more resampling stop helping?
2. Change `confidence` to 0.5 and 0.99. Do the intervals scale the way a normal-theory
   interval would?
3. Fit the three simulated repeats in 4.6 and compare their spread with the bootstrap
   interval. Which is wider, and which one do you believe?
4. Explain, in two sentences, why a bootstrap over 27 rows from nine concentrations is not
   the same as a bootstrap over 27 independent experiments.

## 4.10 What this module did **not** cover

* BCa or studentised bootstrap intervals (DiCiccio & Efron 1996), which correct for bias
  and skewness; the package implements the percentile interval only, and says so.
* Profile-likelihood intervals, which are the best frequentist option for nonlinear models.
* Bayesian posterior intervals, which answer a different question.
* Error propagation for derived quantities (`kcat/Km`): compute it by resampling the
  parameters jointly, do not divide two intervals.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ValueError: the response is constant` | A bootstrap of constant data has no interval to find | Fix the data, not the bootstrap |
| Interval moves every time you run it | You changed the seed without recording it | Report the seed; a reproducible interval is the only kind that counts |
| `ValueError: bootstrapping fewer than four points ...` | Too few points for the refits to be identifiable | More data, or a different method |
| The bootstrap interval is *narrower* than the SE interval | Usually small n and a skewed sampling distribution | Compare both, report the wider one, and say why |

---

**Next:** [Module 5 — weighting, error models and the shape of noise](05-weighting-error-models-and-the-shape-of-noise.md)
