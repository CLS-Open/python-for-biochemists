# Practical 1 — Your first Michaelis–Menten fit

**Time:** 60 minutes. **Needs:** [installation](../installation.md), Module 1.

You will fit Vmax and Km to a committed dataset, check the fit against the simulated truth,
and write the methods sentence you would put in a report.

---

## Step 0 — Look at the file first (5 min)

```bash
head -5 datasets/sample/michaelis_menten_mean.csv
```

```
S,v,sd,n
1.0,0.3601,0.02161,3
2.5,0.92771,0.03747,3
5.0,1.68099,0.01987,3
10.0,2.96149,0.03702,3
```

Answer these before touching any software:

1. What are the units of each column? *(µM, µmol min⁻¹ mg⁻¹, same as v, count — see
   `datasets/metadata/michaelis_menten_mean.yml`.)*
2. How many substrate concentrations? *(9)*
3. What is `sd`? *(sample SD of three replicates; the metadata says so)*
4. Does `sd` grow with `[S]`? *(Not monotonically — check `[S] = 100` against `[S] = 200`
   in the whole file: 0.31285 and 0.25232. The generator mixes a 3 % proportional term with
   an absolute floor, and each `sd` is itself a sample SD of only three replicates, so it is
   a noisy estimate of a quantity that does trend upward. `datasets/metadata/
   michaelis_menten_mean.yml` states the model: `error_model.name:
   proportional_plus_absolute`.)*

## Step 1 — Fit it (10 min)

```bash
kinetics mm datasets/sample/michaelis_menten_mean.csv
```

Record on your worksheet:

| Quantity | Your value |
|---|---|
| Vmax | |
| Vmax standard error | |
| Km | |
| Km standard error | |
| R² | |
| AICc | |
| weighting used | |

Expected (this machine): Vmax **9.76057** (SE 0.168), Km **23.9224** (SE 0.655), method
`weighted_nonlinear`, R² 0.99873, AICc −34.25. The CLI weights by `1/sd²` from the `sd`
column, which is why the estimates are lower than the unweighted ones you will meet in Step 3.
If your numbers differ in the third decimal that is fine; if the *weighting row* differs,
something is wrong with your file.

## Step 2 — Check against the truth (10 min)

The metadata records the parameters used to generate the file:

```bash
grep -A4 true_parameters datasets/metadata/michaelis_menten_mean.yml
```

```
true_parameters:
  vmax: 10.0
  km: 25.0
```

Compute the relative error of each estimate. On this dataset the weighted fit lands at
Vmax −2.39 % and Km −4.31 % from the truth, and the unweighted fit you compute in Step 3 at
−0.16 % and −0.49 %. **You will not have this luxury with real data** — which is exactly why
the exercise exists: it calibrates your sense of what "close" looks like when the error model
is known.

## Step 3 — Same data, no weighting (10 min)

```bash
python - <<'PY'
import csv
from kinetics.fitting import fit_michaelis_menten
rows = list(csv.DictReader(open("datasets/sample/michaelis_menten_mean.csv")))
S = [float(r["S"]) for r in rows]; v = [float(r["v"]) for r in rows]
sd = [float(r["sd"]) for r in rows]
for label, kwargs in (("unweighted", {}), ("weighted", {"sigma": sd})):
    f = fit_michaelis_menten(S, v, **kwargs)
    print(label, round(f.params["vmax"], 3), round(f.params["km"], 3), round(f.aicc, 2))
PY
```

Which is closer to the truth of (10.0, 25.0)? *(Unweighted: 9.984 and 24.877 against
9.761 and 23.922 for the weighted fit — the weighted fit's AICc is also 5.88 units worse,
−34.25 against −40.13.)* Write one sentence explaining why weighting is not automatically
better on a single dataset. (Module 5, section 5.2, discusses this exact result.)

## Step 4 — Write the methods sentence (15 min)

Using the four-part template from Module 2 §2.8 (conditions → model → fitting method →
result with uncertainty), write 3–5 sentences describing what you did. Include:

* the unit of rate and the denominator;
* how many concentrations and replicates;
* the model as a formula;
* the weighting and why;
* both estimates with uncertainties.

## Step 5 — Break it on purpose (10 min)

```bash
kinetics mm datasets/sample/malformed_example.csv; echo "exit=$?"
```

1. What exit code did you get, and why? *(2 — the file has a non-numeric cell, and the
   message names the line: `line 4: v='3.4x' is not a number`.)*
2. Fix the file in `/tmp` (replace `3.4x` with `3.4`) and run again. What now?
   *(Exit **0**. It fits. Vmax comes back as 20.84 ± 30.4 and Km as 61.42 ± 132 — a standard
   error larger than the estimate. The only warning is a note in the output: "1 of 6 rates
   are not positive (lowest −0.4)". Nothing refuses the point for you.)*
3. Remove the negative row too and refit. What happens to the answer?
   *(Vmax 9.893 ± 0.19, Km 9.748 ± 0.53 — on a file of the same shape that was quietly
   reporting nonsense one line earlier. The lesson is not "the software caught it"; it is
   that a fit which runs is not a fit that means anything, and the `sd`/SE columns are how
   you notice.)*

## What to hand in

* Your worksheet table from Step 1.
* The relative errors from Step 2.
* Your two-sentence answer from Step 3.
* Your methods paragraph from Step 4.
* Your exit codes and one-line explanation from Step 5.

## Marking guide

| Criterion | Points |
|---|---:|
| File inspected before fitting; units and replicates identified | 3 |
| Fit reproduced; all six quantities recorded | 4 |
| Comparison with simulated truth, with relative errors computed | 3 |
| Weighting compared and explained in one accurate sentence | 4 |
| Methods paragraph contains all five required elements | 5 |
| Exit codes explained correctly (2 = unreadable, 1 = no answer) | 3 |
| **Total** | **22** |

Marks are additive: a wrong number with a right method earns the method marks. A right
number with no stated weighting earns 2 of 4 for that row.

---

**Next:** [Practical 2 — the Lineweaver–Burk trap](practical-02-the-lineweaver-burk-trap.md)
