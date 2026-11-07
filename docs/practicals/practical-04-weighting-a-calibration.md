# Practical 4 — Weighting a calibration, and an LOD that moves by two orders of magnitude

**Time:** 75 minutes. **Needs:** Module 5.

A detection limit is quoted to three significant figures in papers, certificates and
regulatory filings. In this session you will compute one from the same seven standards four
different ways and watch it move from **0.109** to **11.7** — a factor of 107 — without any
number in the underlying data changing.

---

## Step 1 — Two files, same standards (10 min)

```bash
head -4 datasets/sample/calibration_standards.csv
head -4 datasets/sample/calibration_mean.csv
```

```
concentration,signal,sd
0.0,0.01242,0.001641
5.0,0.245111,0.007319
```

`calibration_standards.csv` has 21 rows (7 concentrations × 3 replicates) and **no sd
column**. `calibration_mean.csv` has 7 rows and carries the replicate SD per standard.

Three questions before you run anything:

1. Which file has more information? *(The 21-row file — unless you know how the replicates
   were drawn, in which case the 7-row file carries almost all of it plus an estimate of the
   noise at each level.)*
2. What does the `sd` column let a fit do that the raw rows do not?
3. Why is the SD at `concentration = 0` (0.001641) the smallest in the file?
   *(It is the blank — the detector floor with no analyte present. Keep this number in mind;
   it is the one that decides the answer in Step 3.)*

## Step 2 — The unweighted calibration (10 min)

```bash
kinetics calibrate datasets/sample/calibration_standards.csv
```

Record: slope **0.0476353**, intercept **0.0124528**, slope SE 0.0005475, residual SD
**0.168369**, n = 21, weighted **False**, LOD **11.663978**, LOQ **35.345387**.

The tool reports `weighted False` because there is no `sd` column to weight by. So it
assumed constant variance — one σ for every standard — and estimated that σ as the residual
SD of all 21 points: 0.168.

That assumption is wrong here, and you can read it straight out of the `sd` column of
`calibration_mean.csv`: 0.001641 at 0 µM, 0.007319 at 5 µM, 0.341712 at 100 µM. The scatter
grows with concentration, because the generator adds 4 % proportional error on top of a
0.0025 blank floor. Ordinary least squares gives the high-concentration points, which carry
that large absolute scatter, exactly the same influence as the blank — which is the wrong
way round when the low end is what you are trying to characterise.

## Step 3 — The same standards, weighted (15 min)

```bash
kinetics calibrate datasets/sample/calibration_mean.csv
```

Record: slope **0.0490126**, intercept **0.0123425**, slope SE **0.0003574**, residual SD
0.0016233, n = 7, weighted **True**, LOD **0.109296**, LOQ **0.331201**.

1. The LOD fell from 11.664 to 0.109 — a factor of **106.7** — while the slope moved by 1.0 %.
   Where did the factor come from? Read the tool's notes: in the weighted fit, σ is
   **0.001641**, the SD at the lowest standard, measured in signal units. In the unweighted
   fit, σ was the residual SD of a fit dominated by the high standards, 0.168369. The ratio
   of those two σ values — 102.6 — is essentially the whole factor, because
   `LOD = 3.3σ/S` is linear in σ.
2. The residual SD fell from 0.168369 to 0.0016233. How much of that is just averaging three
   replicates? *(√3 = 1.732; the observed ratio between the 21-row unweighted residual SD and
   the 7-row unweighted residual SD is 0.168369 / 0.091774 = 1.83. Almost all of it. This is
   why the two numbers are not comparable as "noise levels" — one is scatter about the line,
   the other is scatter of means about the line.)*

## Step 4 — Isolate weighting from n (15 min)

The comparison in Step 3 changed two things at once: the σ definition **and** the number of
rows. Fix that. Strip the `sd` column and refit the seven means unweighted:

```bash
cut -d, -f1,2 datasets/sample/calibration_mean.csv > /tmp/cal_unweighted.csv
kinetics calibrate /tmp/cal_unweighted.csv
```

Record: slope **0.0476353**, residual SD **0.0917743**, n = 7, weighted **False**, LOD
**6.357793**, LOQ **19.266039**.

Now the only difference between this and Step 3 is the weighting:

| Fit | Rows | σ used | Slope | LOD |
|---|---:|---|---:|---:|
| unweighted, means | 7 | residual SD 0.0917743 (signal units) | 0.0476353 | 6.357793 |
| weighted, means | 7 | SD at lowest standard 0.001641 (signal units) | 0.0490126 | 0.109296 |

**LOD ratio: 58.17.** This is the number quoted on the repository home page and in Figure 4.
The simulated truth for the slope is **0.0485**: the weighted estimate is **+1.06 %** away,
the unweighted one **−1.78 %** away.

Answer:

1. Why is the correct σ for an LOD at zero concentration *not* the residual SD of a fit
   spanning 0–200 µM?
2. The weighted fit's notes mention a "scale factor" of 0.989. What is it, and what would a
   value of, say, 4.0 have told you? *(It compares the weights you supplied with the scatter
   actually observed: 1.0 means the declared SDs matched the residuals. A factor of 4 means
   your declared σ were half the true scatter, which usually means the SD column was computed
   from too few replicates or after some averaging.)*
3. Which LOD would you put in a paper, and what exactly must you state alongside it?

## Step 5 — Write the paragraph (25 min)

Write 4–6 sentences for a methods section reporting this calibration and its LOD. Include:

* the calibration range, the number of standards and replicates;
* the model (`signal = a + b·[analyte]`) and the fitting method (weighted least squares by
  `1/sd²`);
* **the σ you used for the LOD, with its definition** — "the standard deviation of the blank"
  is not the same claim as "the residual SD of the regression", and a reader cannot reproduce
  your number without knowing which you meant;
* the LOD with its multiplier (`3.3σ/S`) and the LOQ (`10σ/S`);
* one sentence on what the LOD is *not*: a concentration you can quantify. Three times the
  blank SD is where you can just about tell a signal from noise; the LOQ at **0.331201** is
  where the relative uncertainty is fit for a measurement.

## What to hand in

* The recorded outputs of Steps 2, 3 and 4.
* Your four-answer set from Step 3 and Step 4.
* The completed comparison table.
* The Step 5 paragraph.

## Marking guide

| Criterion | Points |
|---|---:|
| Step 2 and 3 outputs recorded with the weighted flag noticed | 4 |
| Step 3 explains the factor of ~107 via the two different σ values | 6 |
| Step 3 answer 2 identifies replicate averaging (√3) | 2 |
| Step 4 table complete and the LOD ratio computed | 4 |
| Step 4 answer 2 explains the scale factor | 2 |
| Step 5: σ definition stated, LOD and LOQ distinguished, range and n given | 7 |
| **Total** | **25** |

Quoting "LOD = 0.11 µM" without saying which σ produced it earns 2 of 7 for Step 5. The
number is not the result; the number *and its definition* are the result.
