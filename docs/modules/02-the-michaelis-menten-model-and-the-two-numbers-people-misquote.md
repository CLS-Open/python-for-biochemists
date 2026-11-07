# Module 2 — The Michaelis–Menten model and the two numbers people misquote

**Time:** 90 minutes. **Prerequisite:** [Module 1](01-units-rates-and-what-an-assay-measures.md).

## Learning objectives

By the end of this module you can:

1. write the Michaelis–Menten equation and state the three assumptions behind it;
2. explain why Vmax is an extrapolation and compute the substrate concentration needed to
   measure any given fraction of it;
3. fit Vmax and Km and read the two checks the fit needs — the standard errors and the residual
   pattern — rather than R² alone;
4. recognise a fit that is precise and wrong, using the Hill dataset as the demonstration.

---

## 2.1 The equation, and the assumption hidden in it

```
v = Vmax · [S] / (Km + [S])
```

The derivation most courses give (Michaelis & Menten 1913; the Briggs–Haldane
steady-state version) requires:

1. a single substrate, single product reaction;
2. enzyme concentration far below substrate (`[E]₀ ≪ [S]₀`);
3. the steady-state approximation: `d[ES]/dt ≈ 0` after a short pre-steady-state;
4. initial rates, so that reverse reaction and product inhibition are negligible;
5. a single enzyme species with one binding site (or independent, equivalent sites).

Every one of those can fail in a real assay. When they fail, the fitted "Km" still comes
out with a standard error — the fit does not know that its model is wrong. That gap between
*precision* and *correctness* is the theme of the rest of the course.

## 2.2 What Km is

`Km` is the substrate concentration at which the rate is half of Vmax. That is a
definition, and it is the only one that survives contact with data: `Km = (k₋₁ + k₂)/k₁`
under the Briggs–Haldane assumptions, so Km equals the dissociation constant `Kd` only when
catalysis is much slower than substrate release (`k₂ ≪ k₋₁`). For many enzymes it is not.

Consequences worth memorising:

* **Km is an operational constant**, measured under stated conditions. It is not
  automatically a binding affinity, and quoting it as one is a mechanism claim smuggled
  into a number.
* **Km has units of concentration** and does not depend on how much enzyme was used. That
  is why it travels between papers better than Vmax — and why comparisons across buffers,
  pH, or temperatures are so often the source of "disagreement" in the literature.
* **Km is not a threshold.** At `[S] = Km` the rate is half of Vmax; at `Km/10` it is 9 %
  of Vmax; at `10·Km` it is 91 %.

## 2.3 What Vmax is not

Vmax is the limit of the rate as `[S] → ∞`. No experiment reaches it. Every "Vmax" you have
ever seen is an extrapolation from finite data, and its uncertainty is dominated by the
highest substrate concentration you actually measured.

The library says this in the result itself, not only in a manual:

```python
import csv
from kinetics.fitting import fit_michaelis_menten

rows = list(csv.DictReader(open("datasets/sample/michaelis_menten_mean.csv", encoding="utf-8")))
S = [float(r["S"]) for r in rows]
v = [float(r["v"]) for r in rows]
fit = fit_michaelis_menten(S, v)
print(f"Vmax = {fit.params['vmax']:.3f} +/- {fit.stderr['vmax']:.3f}")
print(f"Km   = {fit.params['km']:.3f} +/- {fit.stderr['km']:.3f}")
print("R2 =", round(fit.r_squared, 5), " AICc =", round(fit.aicc, 2))
print(fit.notes[0])
```

```
Vmax = 9.984 +/- 0.102
Km   = 24.877 +/- 0.802
R2 = 0.99934  AICc = -40.13
Vmax is an extrapolation to infinite [S] and is never measured directly
```

```
Vmax is an extrapolation to infinite [S] and is never measured directly
```

(Your digits may differ in the third decimal on a different SciPy build; the note will not
change.)

`kcat` — the turnover number — is `Vmax` divided by the concentration of active sites, so
it inherits the extrapolation *and* adds the active-site titration problem. See
[limitations](../limitations.md) for why this repository does not compute `kcat` for you.

## 2.4 The starting-point problem nobody mentions

Nonlinear least squares finds a *local* minimum. What you get depends on where the
optimiser started, and on how the residuals are weighted. The package makes the starting
point explicit and documented rather than tuned per dataset:

```
the default start is (max observed rate, substrate concentration nearest half of it)
```

That is crude, monotone, and reproducible. If it fails to converge, the fit raises
(`ValueError: fit did not converge for michaelis_menten: ...`) instead of quietly retrying
with different starts, because a fit that only succeeds from a lucky starting point is a
fragile method that should be reported as such.

You can pass your own start:

```python
fit = fit_michaelis_menten(S, v, p0=(12.0, 40.0))
```

## 2.5 Your first fit, and the two checks it needs

```bash
kinetics mm datasets/sample/michaelis_menten_mean.csv
```

The table gives Vmax = 9.98 with a standard error, Km = 24.88, R² = 0.99873 and
AICc = −34.25 for the weighted fit. Two checks before believing any of it:

1. **Do the residuals have a shape?** A curve that sits above the data at low `[S]` and
   below it at high `[S]` means the model is wrong, whatever R² says.
2. **Is the weighting defensible?** This run used `1/sd²` because the file has an `sd`
   column. Without one it would have been unweighted, and Module 5 explains why that is a
   decision rather than a default.

## 2.6 A fit that is precise and wrong: fitting Michaelis–Menten to sigmoid data

`datasets/sample/cooperative_hill.csv` is a sigmoid series generated with `n = 2`
(recorded in its metadata). Fit Michaelis–Menten to it anyway:

```python
import csv
from kinetics.fitting import fit_michaelis_menten, fit_nonlinear
from kinetics.models import hill

rows = list(csv.DictReader(open("datasets/sample/cooperative_hill.csv", encoding="utf-8")))
S = [float(r["S"]) for r in rows]
v = [float(r["v"]) for r in rows]

mm = fit_michaelis_menten(S, v)
h = fit_nonlinear(S, v, hill, (12.0, 30.0, 1.0), ("vmax", "k_half", "n"), name="hill")
print(
    f"Michaelis-Menten: Vmax {mm.params['vmax']:.2f}, Km {mm.params['km']:.2f}, "
    f"R2 {mm.r_squared:.4f}, AICc {mm.aicc:.2f}"
)
print(
    f"Hill:             Vmax {h.params['vmax']:.2f}, K1/2 {h.params['k_half']:.2f}, "
    f"n {h.params['n']:.2f}, R2 {h.r_squared:.4f}, AICc {h.aicc:.2f}"
)
```

```
Michaelis-Menten: Vmax 15.80, Km 53.22, R2 0.9502, AICc 4.37
Hill:             Vmax 11.66, K1/2 28.97, n 2.04, R2 0.9978, AICc -22.59
```

The measured truth is `Vmax = 12.0`, `K½ = 30.0`, `n = 2.0`. Forcing a hyperbola onto
sigmoid data inflates Vmax by 32 % and Km by 77 % — and R² is still 0.95, which is why R²
is a terrible model-selection criterion. Module 6 turns this comparison into a rule
(ΔAICc = 26.97 here).

## 2.7 When the fit refuses

The package refuses in four situations, each with a different message, and each is worth
recognising:

| Situation | What you see | What it means |
|---|---|---|
| Non-numeric cell in the file | `line 3: v='3.4x' is not a number`, exit 2 | fix the file |
| Negative rate | `negative rates cannot be fitted by these mechanisms` | the assay or the blank subtraction is wrong |
| Fewer points than parameters | `need more points (2) than parameters (2)` | collect data; this is not a software problem |
| No convergence | `fit did not converge for michaelis_menten: ...` | the data may not identify the model |

## 2.8 What to write in a methods section

A complete sentence, in this shape:

> Initial rates were measured at nine substrate concentrations (1–200 µM, three replicates)
> at 25 °C in 50 mM Tris-HCl pH 7.4, with total protein 0.4 mg mL⁻¹. Rates
> (µmol min⁻¹ mg⁻¹) were fitted to `v = Vmax·[S]/(Km + [S])` by weighted nonlinear least
> squares, weighting each point by the inverse variance of its replicates. The fit gave
> Vmax = 9.98 ± 0.16 and Km = 24.9 ± 1.2 µM (asymptotic standard errors, n = 9); a
> 95 % bootstrap interval for Km (2000 resamples, seed 20261107) was 23.4–26.4 µM.

Everything in that sentence is a decision you can now defend, including the seed.

## 2.9 Exercises

1. Compute the fraction of Vmax reached at `[S] = Km/10`, `Km`, `10·Km`.
2. Explain why Vmax's standard error is dominated by the highest `[S]` you measured.
3. `datasets/sample/michaelis_menten.csv` has `[S]` up to 200 µM with Km ≈ 25 µM. If you
   truncated the series at 60 µM, would Vmax or Km be hurt more? Test it.
4. The Hill fit reports `n = 2.04 ± 0.11`. Does that prove cooperativity?

   *(No. It shows the curve is steeper than a hyperbola. A monomer with a bad assay can
   produce the same shape; Module 6 explains what would be needed to claim a mechanism.)*

## 2.10 What this module did **not** cover

* Mechanisms with two substrates (Cleland notation), product inhibition, or dead-end
  complexes.
* Allosteric models with physical meaning (MWC, KNF): the Hill equation is a shape
  description, not a mechanism.
* Tight-binding inhibition.
* Time-dependent inactivation.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Km comes back larger than the highest [S] | The substrate range does not bracket Km, so it is extrapolated | Extend the range, or report that Km is not determined |
| Every standard error is `nan` | SciPy could not estimate the covariance (over-parameterised or degenerate data) | Read the note in the output; reduce the model or add data |
| Fit returns values far from a sensible guess | The starting guess `p0` came from data with no curve in it | Supply `p0` explicitly and say so; a converged fit from a lucky start is not a result |
| `ValueError: need more points (5) than parameters (2)` | Fewer rows than parameters | Get more data; the message is the whole story |

---

**Next:** [Module 3 — linearisations and what the transformation costs](03-linearisations-and-what-the-transformation-costs.md)
