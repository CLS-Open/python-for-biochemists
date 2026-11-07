# Practical 3 — Bootstrap intervals, and the difference between replicates and evidence

**Time:** 60 minutes. **Needs:** Module 4.

You will compute a confidence interval three ways — normal theory, resampling the 27
replicates, resampling the 9 means — and discover that the three agree less than the
textbook suggests. Then you will write the sentence that keeps you honest about which one
you are entitled to report.

---

## Step 1 — One fit, two kinds of number (10 min)

```bash
kinetics mm datasets/sample/michaelis_menten.csv
```

`datasets/sample/michaelis_menten.csv` is the full series: 9 substrate concentrations,
3 replicates each, 27 rows. Record:

| Quantity | Value |
|---|---|
| n | |
| Vmax ± SE | |
| Km ± SE | |
| AICc | |

Expected (this machine): n = 27, Vmax **9.9838 ± 0.1069**, Km **24.8774 ± 0.84**, AICc
−97.06.

Now compute the normal-theory 95 % interval by hand, `estimate ± 1.96·SE`. That is the
interval most software prints when it says "95 % CI".

## Step 2 — Bootstrap the same fit (10 min)

```python
import csv
from kinetics.fitting import fit_michaelis_menten, bootstrap_ci

rows = list(csv.DictReader(open("datasets/sample/michaelis_menten.csv")))
S = [float(r["S"]) for r in rows]
v = [float(r["v"]) for r in rows]

ci = bootstrap_ci(S, v, n_resamples=2000, seed=20261107)
print(ci)
```

Expected: Vmax **9.7095 – 10.2630**, Km **23.3591 – 26.4488**.

Compare with `9.9838 ± 1.96 × 0.1069` = 9.774 – 10.193, a half-width of 0.2095 on each side.
The bootstrap interval is **wider at both ends** — 0.2743 below the estimate and 0.2792 above —
and only very slightly asymmetric. The size of the widening is the point: the normal-theory
interval is about 31 % too narrow here, because Vmax is a ratio and its sampling distribution
has a longer upper tail than a bell. On thinner data the asymmetry grows; this dataset is good
enough that it barely shows.

## Step 3 — Nine means, not 27 pieces of evidence (15 min)

The 27 rows are not 27 independent experiments. They are 9 conditions measured three times.
Fit the means instead, and bootstrap them:

```python
means = list(csv.DictReader(open("datasets/sample/michaelis_menten_mean.csv")))
Sm = [float(r["S"]) for r in means]
vm = [float(r["v"]) for r in means]

print(fit_michaelis_menten(Sm, vm).params)
print(bootstrap_ci(Sm, vm, n_resamples=2000, seed=20261107))
```

Expected: the point estimates are **identical to four decimal places** (Vmax 9.9838,
Km 24.8774 — each mean is the same three numbers averaged, so the fit lands in the same
place), but the interval is **wider**: Vmax **9.6349 – 10.3246**, width 0.690 against 0.554
for the 27-row version.

Answer these:

1. Tripling the number of rows narrowed the interval from 0.690 to 0.554 — a 20 % reduction.
   Under the independence assumption, going from 9 to 27 observations should narrow an
   interval by `1 − 1/√3` = 42 %. Why didn't it?
   *(The replicates are repeats of the same 9 conditions. They estimate the within-condition
   noise very well, but they add almost no new information about where the curve bends. An
   interval computed from rows that treats them as independent is an interval about
   *replicate* variability, not about *condition* variability — the pseudoreplication that
   Hurlbert wrote about, listed in [the references](../references.md).)*
2. Which of the two intervals would you report for "how well do we know Km here", and which
   would you report for "how variable is a single assay"?
3. The AICc values are −97.06 for 27 rows and −40.13 for 9 rows. Can you compare them?
   *(No. AICc is comparable only between models fitted to the same response data. Changing n
   changes every term in the formula. The `report` module prints "not comparable" rather
   than a number when the responses differ — that behaviour exists because this comparison
   was made by mistake during development.)*

## Step 4 — Is the bootstrap reproducible? (10 min)

Run Step 2 with five different seeds (`seed=1` … `seed=5`) and record the Vmax interval each
time.

| Seed | Vmax interval | Width |
|---:|---|---:|
| 1 | | |
| 2 | | |
| 3 | | |
| 4 | | |
| 5 | | |

Expected (this machine): 9.709–10.265, 9.708–10.239, 9.700–10.265, 9.710–10.270,
9.715–10.252; widths 0.556, 0.531, 0.565, 0.559, 0.536.

The endpoints move by about 0.03 units between seeds. That is the Monte-Carlo error of the
method, not a property of your data — which is exactly why `bootstrap_ci` takes a **required**
`seed` argument and records it in the result. An interval you cannot reproduce is not a
result you can put in a paper.

## Step 5 — What the bootstrap cannot fix (15 min)

Write 4–6 sentences answering a colleague who says: *"I bootstrapped it, so the error bars
are assumption-free."*

Your answer must include:

* one number from your own Step 2 or 3 output;
* the fact that resampling rows cannot create information about a systematic error shared by
  all replicates (a mis-calibrated pipette resamples perfectly with itself, so its effect
  vanishes from the interval);
* the assumption the bootstrap **does** still make (that the observed points are a
  representative sample of the population the question is about);
* one thing it does not give you at all: any warning about model misspecification. If the
  wrong model is fitted, the bootstrap will produce a tight interval around the wrong answer.

Then state which interval from Steps 1–3 you would put in a methods section, in one sentence,
with the reason.

## What to hand in

* The Step 1 table and your hand-computed normal-theory interval.
* The Step 2 and Step 3 intervals, with the width comparison.
* Your answers to Step 3's three questions.
* The Step 4 table.
* The Step 5 paragraph and your one-sentence methods choice.

## Marking guide

| Criterion | Points |
|---|---:|
| Step 1 table correct and normal-theory interval computed | 3 |
| Step 2 bootstrap reproduced and asymmetry noticed | 4 |
| Step 3 identifies pseudoreplication rather than blaming randomness | 5 |
| Step 3 answer 3: says the AICc values are not comparable | 2 |
| Step 4 table complete, and the seed-dependence is explained | 2 |
| Step 5: shared-error limitation, model-misspecification limitation, defensible choice | 4 |
| **Total** | **20** |

A Step 5 answer that says "the bootstrap is assumption-free" without the shared-error point
scores 1 of 4, however confident it sounds. The interval is assumption-light, not
assumption-free, and the difference is the whole lesson.
