# Practical 2 — The Lineweaver–Burk trap, reproduced by you

**Time:** 75 minutes. **Needs:** Module 3.

You have been told the double-reciprocal plot is biased. This practical makes you produce
the bias yourself, so that the claim is yours rather than the documentation's.

---

## Step 1 — One dataset, four answers (10 min)

```bash
kinetics compare datasets/sample/michaelis_menten_mean.csv
```

From the parameters table, record Vmax and Km for each method, and compute the relative
error against the truth (10.0, 25.0 µM):

| Method | Vmax | rel. error | Km | rel. error |
|---|---:|---:|---:|---:|
| michaelis_menten | | | | |
| hanes | | | | |
| eadie_hofstee | | | | |
| lineweaver_burk | | | | |

Expected: Lineweaver–Burk gives Vmax **11.9869** (+19.9 %) and Km **31.9047** (+27.6 %) on
the same nine points where the nonlinear fit returns 9.98381 (−0.16 %) and 24.8774
(−0.49 %). Hanes–Woolf and Eadie–Hofstee land within 0.7 % on both parameters — the damage
is specific to the double reciprocal, not to linearisation in general.

## Step 2 — Why it happened: transform one point by hand (15 min)

At `[S] = 1 µM` the mean rate is 0.3601. Its replicate SD is 0.02161.

1. Compute `1/v` and `1/S`.
2. Propagate the SD: for `y = 1/v`, `σ_y ≈ σ_v / v²`. Compute `σ_v/v²`.
3. Compare with the transformed value. What fraction of `1/v` is its own error?
4. Repeat for `[S] = 200 µM` (v = 8.86848, SD = **0.25232** — read the row yourself;
   0.31285 belongs to `[S] = 100 µM`).

*Check your answer: `1/v` is 2.77701 with σ 0.16665 at 1 µM, a 6.0 % error fraction; at
200 µM it is 0.11276 with σ 0.00321, a 2.8 % error fraction. The transformed value at the
lowest concentration carries more than twice the relative error of the one at the highest.*

Write down the two error fractions. This is the whole mechanism of the bias: the
double-reciprocal plot weights the point at 1 µM far more heavily than the point at 200 µM,
although the 1 µM point's *relative* error is the same and its absolute rate is 25× smaller.

## Step 3 — Reproduce the distribution (25 min)

One dataset can be lucky. Run the experiment over 200 simulated datasets:

```python
import numpy as np
from kinetics.fitting import compare_linearisations
from kinetics.models import michaelis_menten

S = np.array([1.0, 2.5, 5.0, 10.0, 20.0, 40.0, 60.0, 100.0, 200.0])
labels = ["michaelis_menten", "hanes", "eadie_hofstee", "lineweaver_burk"]
records = {label: [] for label in labels}
for seed in range(200):
    rng = np.random.default_rng(seed)
    # the error model the committed datasets use (see datasets/metadata/*.yml):
    # 3 % proportional error plus an absolute detector floor of 0.05
    v = michaelis_menten(S, 10.0, 25.0) * (
        1 + 0.03 * rng.standard_normal(S.size)
    ) + 0.05 * rng.standard_normal(S.size)
    if np.any(v <= 0):
        continue
    for label, result in zip(labels, compare_linearisations(S, v), strict=True):
        records[label].append(result.params["vmax"])

for label in labels:
    values = np.array(records[label])
    print(
        f"{label:16} median {np.median(values):6.3f}"
        f"  IQR {np.percentile(values, 75) - np.percentile(values, 25):6.3f}"
        f"  within 20% {np.mean(np.abs(values - 10) < 2):5.0%}"
    )
```

Fill in:

| Method | median Vmax | IQR | fraction within 20 % |
|---|---:|---:|---:|
| nonlinear | | | |
| Hanes–Woolf | | | |
| Eadie–Hofstee | | | |
| Lineweaver–Burk | | | |

Expected from those 200 datasets: nonlinear median 9.994 with IQR 0.378, Hanes–Woolf 9.974
with 0.430, Eadie–Hofstee 9.924 with 0.460, and Lineweaver–Burk 9.669 with IQR **4.405**,
only 50 % of estimates within 20 % of the truth, and a worst case of 44.47. Those are the
numbers Figure 2 in the [repository home page](../index.md) is drawn from: nonlinear and
Hanes tight, Eadie–Hofstee looser, Lineweaver–Burk with a long tail.

## Step 4 — Find a dataset where Lineweaver–Burk looks perfect (10 min)

Search seeds until you find one where `lineweaver_burk` recovers Vmax within 1 % of 10.0.

1. How many seeds did you have to try?
2. Print that dataset's Lineweaver–Burk Vmax alongside the nonlinear fit.
3. Write two sentences explaining why "my Lineweaver–Burk plot looked fine" is not
   evidence, using the number you just found.

*Check your answer: across these 200 datasets only 5 recover Vmax within 1 % by
Lineweaver–Burk (the first is seed 62, giving 9.945). A plot that "looks fine" therefore
happens about 2.5 % of the time by luck, which is not a standard anyone can report.*

## Step 5 — What would you tell your supervisor? (15 min)

Write 4–6 sentences for a colleague who plans to publish Km values from double-reciprocal
plots. You must include:

* the mechanism (error amplification by `1/v²`);
* one number from your own Step 3 output;
* a concrete alternative (`fit_michaelis_menten`);
* one legitimate use of the plot (a diagnostic, e.g. Eadie–Hofstee for detecting curvature).

## What to hand in

The four tables, your two error fractions from Step 2, and the paragraph from Step 5.

## Marking guide

| Criterion | Points |
|---|---:|
| Step 1 table complete with relative errors | 4 |
| Step 2 propagation computed correctly for both points | 5 |
| Step 3 script runs and table filled with measured values | 6 |
| Step 4 dataset found and explained | 3 |
| Step 5 paragraph: mechanism + own number + alternative + legitimate use | 7 |
| **Total** | **25** |

A paragraph that says "Lineweaver–Burk is bad" without a mechanism earns 2 of 7, however
strongly worded it is.
