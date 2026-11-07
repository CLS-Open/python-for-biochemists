# Practical 5 — Typing an inhibitor, or refusing to

**Time:** 90 minutes. **Needs:** Module 6.

"What kind of inhibitor is it?" is the question every reviewer asks and almost nobody
answers properly. You will fit the three classical mechanisms to two datasets — one designed
well enough to answer the question, one not — and write the sentence that says which
situation you are in.

---

## Step 1 — The well-designed experiment (15 min)

```bash
kinetics inhibition datasets/sample/inhibition_series.csv
```

Record the model comparison table:

| Model | R² | AICc | ΔAICc |
|---|---:|---:|---:|
| competitive | | | |
| noncompetitive | | | |
| uncompetitive | | | |

Expected (this machine): competitive **−83.17** (Δ 0.00), noncompetitive **−40.62**
(Δ 42.55), uncompetitive **−25.27** (Δ 57.90).

And the parameters of the winner:

```
competitive  vmax 9.6056 ± 0.145   km 24.0369 ± 1.01   ki 8.59114 ± 0.549
```

The simulated truth is Vmax 10.0, Km 25.0 µM, **Ki 8.0 µM, competitive**.

Answer:

1. Which parameter moved furthest from truth, and by how much?
   *(Vmax, by 3.9 %: 9.6056 against 10.0. Km moved 3.9 % the other way, and Ki is 7.4 % high.
   Competitive inhibition raises the apparent Km and leaves Vmax unchanged — but the fit does
   not know that, and with 20 points it trades a little Vmax against a little Km.)*
2. ΔAICc of 42.55 for the next-best model. Is that close? *(Δ > 10 is conventionally treated
   as no support at all. There is no ambiguity in this dataset.)*
3. The tool prints the caveat "a best-fitting mechanism is not a proven mechanism" in the
   header. Why, given the fit is this clear? *(Because the data were generated from a
   competitive model. Real data can be fitted by a mechanism that is not the one operating —
   and the fit cannot see a mechanism it was not given.)*

## Step 2 — The badly designed experiment (20 min)

```bash
kinetics inhibition datasets/sample/inhibition_weak.csv
```

Expected (this machine):

| Model | R² | AICc | ΔAICc |
|---|---:|---:|---:|
| uncompetitive | 0.99755 | −54.09 | 0.00 |
| noncompetitive | 0.99747 | −53.62 | 0.48 |
| competitive | 0.99726 | −52.41 | 1.68 |

1. What is the true mechanism for this dataset? *(Read
   `datasets/metadata/inhibition_weak.yml`: it is **competitive** with **Ki = 500 µM** — a
   weak inhibitor measured with inhibitor concentrations that never go near its Ki.)*
2. Does the best model recover it? *(No. The answer comes back uncompetitive.)*
3. What are the three ΔAICc values, and what do they license you to say?
   *("The three mechanisms are indistinguishable on these data." A Δ of 0.48 between the top
   two is nothing; the third at 1.68 is also inside the "no support for one over another"
   band. The fit is telling you the experiment cannot answer the question.)*
4. Look at the Ki estimates: 47.4696 ± 34.1 for uncompetitive, 51.333 ± 74.9 for competitive,
   against a true Ki of 500 µM. What does a standard error larger than the estimate mean
   here? *(That Ki is not identified. Every Ki in this file is off by an order of magnitude
   and carries an uncertainty wider than itself: the inhibitor series never reached high
   enough [I] to bend the curve, so the fit is extrapolating Ki from noise.)*

## Step 3 — Design the experiment that would answer it (20 min)

`datasets/metadata/inhibition_weak.yml` records how the file was generated, including the
inhibitor concentrations used. Read it, then write down the design you would have used
instead. Your answer must state:

* how many inhibitor concentrations, and at what values relative to Ki (a rule of thumb:
  span 0 to ~10 × Ki, with at least four non-zero levels);
* how many substrate concentrations, and where (at least one well below Km, one at Km, one
  well above);
* the total number of assays, and whether the replicates are of the same cuvette or of the
  whole assay.

Then check your reasoning against the generator:

```bash
grep -A6 "inhibitor_concentrations\|substrate_concentrations" datasets/metadata/inhibition_weak.yml
grep -A6 "inhibitor_concentrations\|substrate_concentrations" datasets/metadata/inhibition_series.yml
```

The contrast is stark and it is in the files:

| Dataset | True Ki | [I] used (µM) | Highest [I] as a fraction of Ki | Outcome |
|---|---:|---|---:|---|
| `inhibition_series` | 8.0 | 0, 2, 8, 16 | 2.0 (200 %) | mechanism typed, ΔAICc 42.55 |
| `inhibition_weak` | 500.0 | 0, 1, 2 | 0.004 (0.4 %) | mechanisms indistinguishable, ΔAICc 1.68 |

Same five substrate concentrations (2, 5, 10, 25, 100 µM) in both files. The only thing that
changed is the inhibitor series — and that single design decision is the difference between a
result and a non-result. `inhibition_weak.yml` says it plainly in `known_problem`: *"[I] never
exceeds 0.004·Ki, so Ki is not identified"*.

## Step 4 — The global fit trap (15 min)

The four-parameter mixed model (Vmax, Km, Ki, **α**) is what most fitting software offers as
"mixed inhibition". This package deliberately does **not** include it in
`compare_inhibition_models`, but it gives you the forward model and the responsibility. Fit it
yourself:

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

Expected output on the committed data:

```
vmax     9.606 ±     nan
km      24.037 ±     nan
ki       8.591 ±     nan
alpha  97030109.768 ±     nan
AICc -80.01
```

Answer:

1. Why is the mixed model absent from `compare_inhibition_models`, and what does that cost
   you? *(Competitive, uncompetitive and noncompetitive are all special cases of it (α → ∞,
   Ki' → ∞, α = 1), so on any dataset it can fit at least as well — and on this design the
   extra parameter costs more AICc than it buys: −80.01 against the competitive fit's
   −83.17. Offering it in the table would invite a number with no interval to be reported as
   a result.)*
2. Read the α estimate: 9.7 × 10⁷ with a `nan` standard error. What has the optimiser
   actually told you? *(α → ∞ is the competitive limit: with four parameters it rediscovered
   the three-parameter answer. And `nan` for every standard error means the covariance matrix
   was not estimable at all — an estimate without an interval is not a measurement.)*
3. What would you have to see to justify reporting α? *(A finite α whose interval excludes
   both limits α → 0 and α → ∞. On a design with [I] on both sides of Ki, that is achievable;
   on this one it is not.)*

## Step 5 — The paragraph (20 min)

Write 5–7 sentences for a paper reporting both datasets. Include:

* the datasets and their designs (numbers of [S] and [I] levels);
* for the well-designed dataset: the winning mechanism, ΔAICc to the next best, and Ki with
  its uncertainty;
* for the weak dataset: that the mechanisms were indistinguishable, the ΔAICc values, and
  explicitly that the fit's preferred mechanism is **not** being claimed;
* one sentence on what experiment would settle it.

## What to hand in

* The two comparison tables.
* Your answers to Steps 1, 2 and 4.
* Your Step 3 design, with the numbers.
* The Step 5 paragraph.

## Marking guide

| Criterion | Points |
|---|---:|
| Step 1 table complete and ΔAICc interpreted correctly | 4 |
| Step 2: recognises the wrong mechanism is reported and says why | 6 |
| Step 2 answer 4: Ki is not identified | 3 |
| Step 3: a design with justified [I] and [S] choices and a stated total n | 6 |
| Step 4: explains the omission of the mixed model, reads α, and says what would justify it | 4 |
| Step 5: both datasets reported, including the explicit non-claim | 7 |
| **Total** | **30** |

A Step 5 paragraph that reports "uncompetitive inhibition" for `inhibition_weak.csv` without
the non-claim scores 1 of 7, even if every number in it is correct. Reporting the number and
the doubt together is the skill this practical exists to teach.
