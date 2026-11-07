# Prerequisites

This is the CLS-Open **intermediate** repository. It assumes the beginner repository
(`python-for-biologists`) or an equivalent course, and it assumes some biochemistry.

## What you should already be able to do

| Skill | Why it is needed | If you cannot |
|---|---|---|
| Write and run a Python function with arguments and a return value | every example is one | start with CLS-Open #1 |
| Read a traceback and find the line that raised | the course deliberately shows errors | CLS-Open #1, Module 3 |
| Use `numpy` arrays for simple arithmetic | the models are vectorised | 30 minutes with the numpy absolute-beginners guide |
| Explain what a derivative is | the rate is one | any calculus refresher; you only need the idea of a slope |
| Explain what a standard deviation is | half the course is about uncertainty | you need this one properly, see below |

## Biochemistry you need

* **Michaelis–Menten kinetics**: `v = Vmax·S/(Km + S)`. If you have seen the curve and the
  term `Km`, you have enough to start. Module 2 rebuilds it from the assay upward anyway.
* **What an initial rate is**: the slope of product against time in the first part of the
  reaction, before substrate depletion and product inhibition bend the curve.
* **What an inhibitor does**: reduces the rate, by one of several mechanisms. You do not
  need to know the mechanisms in advance; Module 6 is where they are compared.
* **What a calibration curve is**: standards of known concentration, measured signal, a
  straight line, and a limit of detection. Module 5 and Practical 5 cover this from first
  principles.

## Mathematics you need

You will meet these formulas. You do not have to derive them; you do have to be able to
read them.

| Formula | Where | What it means in words |
|---|---|---|
| `v = Vmax·S/(Km + S)` | Module 2 | rate rises with substrate and flattens at Vmax |
| `1/v = (Km/Vmax)(1/S) + 1/Vmax` | Module 3 | the Lineweaver–Burk transformation: a straight line |
| `S/v = S/Vmax + Km/Vmax` | Module 3 | the Hanes–Woolf transformation |
| `AIC = n·ln(RSS/n) + 2k` | Module 6 | fit quality minus a penalty for parameters |
| `AICc = AIC + 2k(k+1)/(n-k-1)` | Module 6 | the small-sample correction, needed whenever n/k < 40 |
| `LOD = 3.3·σ/slope` | Module 5 | the smallest signal distinguishable from blank noise |

If `AIC` is unfamiliar, that is normal: it is taught here, and Module 6 spends its first
half on what it does and does not mean.

## Statistics you should be honest about

The course's hardest idea is not a formula, it is a distinction:

* a **standard error** estimates how much an estimate would move if the experiment were
  repeated, *assuming the model is correct*;
* a **confidence interval** is a statement about a procedure, not about your parameter;
* a **bootstrap interval** is an estimate of the first, obtained by resampling your own
  data, and it inherits whatever your data do not contain.

You do not need to have studied mathematical statistics. You do need to be willing to stop
saying "the Km is 25" and start saying "the Km is estimated at 24.9 ± 0.8 µM under this
model and this weighting". The whole repository is about that shift.

## Software you need

* Python 3.10 or newer.
* A terminal you can run one command in.
* Any editor. Jupyter is optional; the notebook in `notebooks/` is executed for you.
* No internet connection for the core practicals. Only `scripts/check_references.py`
  needs network access.

## What you do **not** need

* A GPU, a cluster, or a database.
* Linear algebra beyond "a matrix is a grid of numbers".
* Prior experience with `scipy.optimize` — but you will leave with some.
* Any real laboratory data. Every dataset here is simulated and says so in
  `datasets/metadata/`, which is what makes the answers checkable.

## A self-test

You are ready if you can answer these three questions without looking anything up. If you
cannot, the linked pages will fix it.

1. A reaction has `Vmax = 10` and `Km = 25 µM`. What is the rate at 25 µM, and what
   fraction of Vmax is reached at 250 µM?
   *(Answer: 5, and 250/275 = 91 %.)*
2. Your fit reports `Km = 24.9 ± 0.8 µM`. Your colleague reports `Km = 25 µM`. Which of
   you is making a claim the data support?
   *(Answer: the one carrying the uncertainty, provided the model and weighting are stated.)*
3. Why can you never measure `Vmax` directly?
   *(Answer: it is the limit as [S] → ∞; see Module 2.)*

If those felt obvious, start at [Module 1](modules/01-units-rates-and-what-an-assay-measures.md).
If they did not, read [Learning path](learning_path.md) first — it tells you which sections
you can skip and which you cannot.
