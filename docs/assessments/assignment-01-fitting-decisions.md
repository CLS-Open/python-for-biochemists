# Assignment 1 — Fit a curve, then defend it

**Released after:** Module 2. **Due:** after Module 4. **Weight:** 10 % of the module mark.
**Time:** 3–4 hours. **Total: 30 marks.**

You are given a substrate series and a supervisor who wants "the Km" by Friday. This
assignment is what you hand back instead: the number, the interval, and the two sentences that
say what the interval does not cover.

---

## The data

Use `datasets/sample/michaelis_menten.csv` (27 rows: 9 concentrations × 3 replicates) **and**
`datasets/sample/michaelis_menten_mean.csv` (9 means with replicate SD). Do not modify either
file; if you want to remove a point, make a copy and say which point and why.

## Tasks

**1. Fit both files (6 marks).** Run the nonlinear fit on each. Record for both: Vmax, Km,
their standard errors, the method name and the weighting flag the tool reports.

*Expected shape of the answer:* the 27-row and 9-row fits return the same point estimates to
four decimal places (Vmax 9.9838, Km 24.8774), because each mean is the average of the same
three numbers. The standard errors differ slightly (0.1069 against 0.1021 for Vmax).
Explain in two sentences why the point estimates agree and the standard errors do not.

**2. State the uncertainty you would report (6 marks).** Compute the normal-theory 95 %
interval (`estimate ± 1.96·SE`) and the percentile bootstrap interval (2000 resamples, seed
20261107) for the 27-row fit. Report both, and state which you would put in a paper with one
sentence of justification. The bootstrap interval is 9.7095 – 10.2630 for Vmax; the
normal-theory one is 9.774 – 10.193. Marks are for the justification, not the numbers.

**3. Check against the truth (4 marks).** The metadata records Vmax = 10.0, Km = 25.0. Compute
the relative error of each estimate. Then do the same for the Lineweaver–Burk fit
(`kinetics compare` gives 11.9869 and 31.9047) and say which of the four methods was closest to
the truth on this dataset, and whether that ordering is a general result.

**4. What the interval does not cover (8 marks).** Write 4–5 sentences naming three things your
interval does **not** account for, each with a concrete example from this dataset:

* a systematic error shared by all replicates (what happens to the interval if the stock
  concentration were 5 % wrong?);
* replicate structure (what 27 rows mean when there are 9 conditions);
* model choice (what would change if the true mechanism were slightly sigmoid).

**5. The methods paragraph (6 marks).** Write the methods paragraph you would publish,
including: the model as a formula, the fitting method, the weighting and why, n and the
replicate structure, the uncertainty method with its seed, and the units of every quantity.

## Rubric

| Criterion | Marks | Full credit looks like |
|---|---:|---|
| 1 · Both fits recorded accurately, agreement explained | 6 | States that the estimates are identical because averaging is linear, and that the SEs differ because the residual structure differs |
| 2 · Both intervals computed, choice justified | 6 | Chooses the bootstrap and says why (ratio, asymmetric sampling distribution), or chooses the SE interval with an argument for why the assumption holds here |
| 3 · Relative errors computed; ordering interpreted | 4 | Says Hanes–Woolf was closest here *and* that this is one dataset, not a law |
| 4 · Three specific limitations with consequences | 8 | Each limitation is tied to this dataset, not quoted generically |
| 5 · Methods paragraph complete and unambiguous | 6 | A reader could re-run the analysis without asking a question |
| **Total** | **30** | |

## Submission rules

* Every number in your write-up must be traceable to a command or script you submit.
* State the Python version, the `kinetics` version (`kinetics --version`) and any seed you used.
* A submission that quotes a number it cannot regenerate is capped at 15/30, however correct
  the number is.

## Common ways this goes wrong

* Reporting the standard error as "±" without saying whether it is 1 SE or a 95 % interval.
* Quoting Vmax to five significant figures while the interval spans ±0.3.
* Answering task 4 with "the model might be wrong" — that is a category, not a limitation of
  *this* analysis.
