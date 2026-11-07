# Challenge 1 — The dataset that cannot be typed

**Released after:** Module 6. **Weight:** 5 % of the module mark. **Time:** 2–3 hours.
**Total: 20 marks.** Do this alone; discussion of the *approach* is allowed, not of the numbers.

`datasets/sample/inhibition_weak.csv` reports "uncompetitive" as the best mechanism. Its true
mechanism is competitive, with Ki = 500 µM. Both statements are in the repository, so neither is
the challenge.

**The challenge is this: change the experiment so that the answer becomes true, without changing
the enzyme.**

## Tasks

1. **(4 marks)** Reproduce the failure. Run the comparison, record the three AICc values and the
   Δ scale, and state in one sentence what the fit is licensed to claim. Then read
   `datasets/metadata/inhibition_weak.yml` and compute, in numbers, why the design cannot
   answer the question (the ratio of the highest [I] to Ki is the key quantity).

2. **(8 marks)** Design a better experiment and simulate it. Write your own generator (not a
   modification of `scripts/make_datasets.py`) that produces a competitive inhibitor with
   Vmax = 10, Km = 25, Ki = 500 using:
   * the same five substrate concentrations as the committed file;
   * your own inhibitor grid;
   * the documented error model (`proportional_plus_absolute`, cv 0.03, sd_abs 0.05) and a seed
     you report.

   Fit all three mechanisms and show ΔAICc. Your design must recover "competitive" with the true
   mechanism ranked first and ΔAICc > 10 to the runner-up.

3. **(5 marks)** Find the boundary. Repeat the simulation with the highest [I] halved, then
   quartered, and report the ΔAICc and the recovered best model at each step. At what fraction of
   Ki does the design stop working? A table of four rows (including your winning design) is
   enough, but every number must come from your script.

4. **(3 marks)** Write the two sentences a supervisor would need: one saying what the original
   dataset supports, and one saying what your design supports. Neither sentence may contain the
   word "proves".

## Rubric

| Criterion | Marks | Full credit |
|---|---:|---|
| 1 · Failure reproduced and the design diagnosed in numbers | 4 | Names [I]max/Ki = 0.004 and links it to the ΔAICc |
| 2 · Own generator, correct model, mechanism recovered | 8 | Generator committed, seed stated, ΔAICc > 10, Ki within 25 % of 500 |
| 3 · Boundary mapped with real numbers | 5 | A trend across at least three designs, with the failure point identified |
| 4 · Two sentences, neither overclaiming | 3 | "Cannot distinguish" in one, "supports competitive" in the other |
| **Total** | **20** | |

## What will fail

* Editing the committed dataset instead of generating your own. The point is that the *design*
  changes, not the file.
* A design that works only because every rate is near zero at high [I]. Report the lowest
  simulated rate and compare it with the blank SD (0.0025 in the error model); a design whose
  signal has vanished has not solved the problem, it has moved it.
* Reporting "competitive wins" without the ΔAICc that licenses the claim.
