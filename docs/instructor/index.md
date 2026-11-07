# Instructor pack

Everything a teacher needs to run this course, including the parts that usually go wrong.

## Who this is for

Learners who have finished CLS-Open #1 (`python-for-biologists`) or who already write Python
and have met a Michaelis–Menten plot. They do **not** need statistics beyond a mean and a
standard deviation. They do need to be comfortable running a command that prints a table they
must then interpret — the course fails if students read output as an answer key.

Typical cohorts:

| Cohort | Route | Time |
|---|---|---|
| Undergraduates, 3rd year biochemistry | Modules 1–4 + Practicals 1–3 | 8–10 h |
| MSc / graduate rotation | Whole course + Practical 5 + mini-project | 20–24 h |
| Statistics or data-science students | Modules 3–6 + Practicals 2, 3, 4 | 12 h |
| Research staff, one-day workshop | Module 3 + Practical 2 + Module 6 | 5 h |

## Timings, routes and the pre-session checklist

Those live in the [teaching plans](teaching-plans.md): the single-day and two-day formats,
the half-day workshop, cohort routes and the install checklist. Read that page first if you
are delivering the course; this page is about the teaching content.

## The five misconceptions this course is built to break

Each one is addressed by a specific measured result, not by an assertion. When a student
states one of these, send them to the artefact rather than correcting them verbally — the
number is more persuasive than the instructor.

| Misconception | Where it dies |
|---|---|
| "R² = 0.999 means the fit is right" | Practical 2 Step 3: Lineweaver–Burk reaches R² = 0.99875 with Vmax 19.9 % high |
| "More data always narrows the interval" | Practical 3 Step 3: 27 rows give a *narrower* interval than 9 means only by 20 %, not 42 % |
| "The best AICc tells me the mechanism" | Practical 5 Step 2: `inhibition_weak` names the wrong mechanism with Δ = 0.48 |
| "Weighting is a detail" | Practical 4: the same standards give LOD 0.109 or 6.358 depending on σ |
| "The software would have warned me" | Practical 1 Step 5: a negative rate fits happily and returns Vmax 20.84 ± 30.4, with one note and no error |

## What to mark, and what to ignore

Mark the **sentences**, not the numbers. A student who reproduces 9.9838 and writes "the fit is
good because R² is high" has learned nothing; a student who gets 9.97 and writes "the estimate is
0.3 % from truth but the interval spans ±0.2, so the third digit is not real" has learned the
course.

Each practical carries its own additive rubric (20–30 marks). Two cross-cutting criteria are
worth marks in every practical:

* **Units and denominators** stated with every number (1 mark per practical).
* **The non-claim**, wherever the data do not support a conclusion (2 marks in Practicals 3
  and 5).

Suggested weighting for a graded module: Practicals 1–5 at 10 % each, quiz 30 %, mini-project
20 %.

## Common questions from students, with answers that hold up

**"Which fitting method should I use?"** Nonlinear least squares on the untransformed data,
weighted if you have a defensible σ. The linearisations are in the package so you can measure
what they cost, not because they are options.

**"My supervisor says Lineweaver–Burk is fine."** Run Practical 2 Step 3 on their data and show
the distribution, not the single line. The argument is about the tail, and the tail is invisible
in one dataset.

**"Can I just drop the outliers?"** You can, but you must say so and show the fit both ways. The
package will not do it silently, and neither should you.

**"Why is Vmax so unstable?"** Because it is an extrapolation to infinite substrate. Report it,
but report Km — which is measured in the middle of your data — as the better-determined number
when that is what the intervals say.
