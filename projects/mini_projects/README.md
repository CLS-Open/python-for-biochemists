# Mini-projects (4–6 hours each)

Two projects, both of which end in something you can hand to a supervisor. Pick one; do both
if you want the capstone to be easier.

Neither has a worksheet. Each has a **brief**, a **definition of done**, and a warning about the
mistake that is most often made. The numbers you report must come from your own runs, and your
write-up must include the commands that produced them.

---

## Mini-project A — Type an inhibitor you cannot type

**Brief.** `datasets/sample/inhibition_weak.csv` is a real competitive inhibitor (Ki = 500 µM)
measured with inhibitor concentrations that never approach its Ki. `kinetics inhibition` returns
"uncompetitive" as the best model, with ΔAICc = 0.48 to the next.

Your job is to **design and simulate the experiment that would answer the question**, then
demonstrate that it does.

1. Read `datasets/metadata/inhibition_weak.yml` and
   `datasets/metadata/inhibition_series.yml` and state, in numbers, what differs between the
   two designs.
2. Write a generator (a Python script in your own repository, not a modification of this one)
   that simulates a competitive inhibitor with Vmax = 10, Km = 25, Ki = 500, using your chosen
   [S] and [I] grids and the same error model as the committed datasets
   (`proportional_plus_absolute`, cv 0.03, sd_abs 0.05, a fixed seed).
3. Show that your design recovers the mechanism: fit all three models, report ΔAICc for each,
   and state the Ki you recover against the truth of 500.
4. Show that a *smaller* version of your design fails — find the [I] range at which the
   mechanisms become indistinguishable again. This is the part that makes the project useful:
   you are mapping the boundary, not just crossing it.
5. Hand in a figure of ΔAICc against maximum [I]/Ki.

**Definition of done.** A table of ΔAICc for at least four designs, a sentence per design saying
what it establishes, and the Ki from the winning design with its standard error.

**The common mistake.** Choosing [I] up to 10·Ki without checking that the rates stay above the
detector floor; at 500 µM Ki, high [I] makes every rate indistinguishable from zero and the fit
degrades even though the design "looks" better. Report the lowest rate in your simulated data
and compare it with the blank SD.

---

## Mini-project B — Does the LOD you report survive a re-analysis?

**Brief.** Choose the LOD from `datasets/sample/calibration_mean.csv` — 0.109296 weighted,
6.357793 unweighted — and write a short report answering: *which one is the LOD, and what would
it take to change my mind?*

1. Recompute both LODs yourself (`kinetics calibrate`, plus the same computation with the `sd`
   column removed) and confirm the reported numbers to four decimal places.
2. Determine the σ each one used, and say which assumption is right for a colorimetric assay
   with a 4 % proportional error and a 0.0025 blank floor.
3. **Stress the decision.** Perturb the input in three ways and report the LOD each time: (a)
   drop the highest standard; (b) multiply every `sd` by 1.5; (c) refit using only the four
   lowest standards. Present all four LODs in one table.
4. Say which of those perturbations matters most, and why that one is a warning about
   reporting a single number.
5. Finish with the methods paragraph you would actually publish, including the σ definition.

**Definition of done.** A table of five LOD values (baseline plus three perturbations plus one
of your own choosing), the σ used in each case, and a paragraph that a reviewer could not
misread.

**The common mistake.** Reporting the LOD to three significant figures without saying which σ
produced it. The two analyses here differ by a factor of 58; a reader given only "LOD = 0.11 µM"
cannot tell whether they could reproduce it, and usually could not.

---

## Marking (both projects)

| Criterion | Points |
|---|---:|
| Every reported number traceable to a command in the submission | 5 |
| Design or perturbation experiment actually run, with numbers | 8 |
| Correct interpretation of the ΔAICc or LOD evidence, including its limits | 5 |
| A stated limitation of the analysis (not a generic one) | 4 |
| Reproducibility: a fixed seed, a committed generator, the environment stated | 3 |
| **Total** | **25** |

A submission whose numbers cannot be regenerated from its own files scores at most 12,
however good the prose. That rule is the one thing this series will not compromise on.
