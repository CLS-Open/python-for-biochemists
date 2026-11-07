# Challenge 2 — The LOD you would defend under review

**Released after:** Module 5. **Weight:** 5 % of the module mark. **Time:** 2 hours.
**Total: 20 marks.**

A reviewer has written: *"The authors report an LOD of 0.11 µM but do not state how it was
calculated. Please clarify."*

Your job is to write the reply that survives the follow-up question — and to find out, by
computation, how much the answer can move before it stops being defensible.

## Tasks

1. **(4 marks)** Compute the LOD for `datasets/sample/calibration_mean.csv` both ways
   (weighted by the `sd` column, and unweighted by stripping it). Record the slope, the σ, the
   LOD and the LOQ for each. Confirm the factor between the two LODs (it is 58.17 to two
   decimal places) and say which σ produced each.

2. **(8 marks)** Stress the number. Report the LOD for each of these variants, with the σ used:
   * all seven standards, weighted (baseline);
   * the four lowest standards only, weighted;
   * all seven, with every `sd` multiplied by 1.5;
   * all seven, with the highest standard removed.

   You may use `kinetics calibrate` on edited copies of the file (say what you changed) or call
   `kinetics.calibration.fit_calibration` directly. Either way, the four LODs must come from
   code you submit.

3. **(5 marks)** Say which of the four is the one you would report, and defend it against the
   obvious objection to each of the other three. ("Dropping the top standard narrows the
   calibration range; that is a different claim, not a better one.")

4. **(3 marks)** Write the reviewer reply, in 4–6 sentences, containing: the σ definition, the
   range and n, the multiplier used for LOD and LOQ, and one sentence distinguishing the two.

## Rubric

| Criterion | Marks | Full credit |
|---|---:|---|
| 1 · Both LODs with σ definitions | 4 | Names 0.001641 (blank SD) and 0.0917743 (residual SD) explicitly |
| 2 · Four variants computed by code | 8 | Every number traceable; σ reported for each |
| 3 · Choice defended against each alternative | 5 | Engages with the objection rather than dismissing it |
| 4 · Reviewer reply complete | 3 | A reviewer could reproduce the number from the reply alone |
| **Total** | **20** | |

## The trap

It is easy to produce a table of four LODs and then report the smallest, or the one from the
"cleanest" subset. Neither is a decision. The mark is for the argument that connects the σ you
chose to the question the LOD is supposed to answer: *how little analyte can I detect in a
sample like the ones I will actually measure?* That question, not the arithmetic, decides which
of the four numbers belongs in the paper.
