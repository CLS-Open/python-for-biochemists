# Learning path

The course is six modules and five practicals. Read in order the first time: the modules
build a single argument, and Module 6 assumes Modules 3–5 have already made you suspicious
of a tidy R².

## The whole course

| Step | What | Time | You can do afterwards |
|---|---|---:|---|
| [Installation](installation.md) | environment | 0:10 | run `kinetics --version` |
| [Module 1](modules/01-units-rates-and-what-an-assay-measures.md) | units and rates | 1:15 | state what your numbers mean |
| [Practical 1](practicals/practical-01-first-michaelis-menten-fit.md) | first fit | 1:00 | fit Vmax and Km with a stated weighting |
| [Module 2](modules/02-the-michaelis-menten-model-and-the-two-numbers-people-misquote.md) | the model | 1:30 | explain why Vmax is an extrapolation |
| [Module 3](modules/03-linearisations-and-what-the-transformation-costs.md) | linearisations | 1:30 | measure the bias instead of asserting it |
| [Practical 2](practicals/practical-02-the-lineweaver-burk-trap.md) | the trap | 1:15 | show a supervisor the tail |
| [Module 4](modules/04-uncertainty-standard-errors-and-the-bootstrap.md) | uncertainty | 1:30 | report an interval you can defend |
| [Practical 3](practicals/practical-03-bootstrap-intervals-and-replicates.md) | resampling | 1:00 | explain what a bootstrap does not fix |
| [Module 5](modules/05-weighting-error-models-and-the-shape-of-noise.md) | weighting | 1:30 | choose a weight and justify it |
| [Practical 4](practicals/practical-04-weighting-a-calibration.md) | calibration | 1:15 | produce a defensible LOD |
| [Module 6](modules/06-choosing-between-models-without-fooling-yourself.md) | model choice | 1:45 | report "indistinguishable" without embarrassment |
| [Practical 5](practicals/practical-05-typing-an-inhibitor.md) | inhibition | 1:30 | type a mechanism, or say why you cannot |
| [Assessment](assessments/quiz-01-modules-1-6.md) | quiz and rubric | 1:00 | check yourself |
| [Mini-project](https://github.com/cls-open/python-for-biochemists/blob/main/projects/mini_projects/README.md) | 4–6 h | 5:00 | run a complete analysis of your own |

Total: about 20 hours for the guided path, plus the mini-project.

## Shorter routes

### One afternoon (4 hours)

For someone who already fits kinetics curves and wants the honest version.

1. [Module 3](modules/03-linearisations-and-what-the-transformation-costs.md) — why the
   plot you were taught is the weakest of the four methods.
2. [Module 4](modules/04-uncertainty-standard-errors-and-the-bootstrap.md) — what the
   error bars you report actually mean.
3. [Practical 2](practicals/practical-02-the-lineweaver-burk-trap.md) — reproduce the
   bias on your own simulated data.
4. [Module 6](modules/06-choosing-between-models-without-fooling-yourself.md) — AICc, and
   when to refuse to choose.

### The statistics route

Skip the biochemistry, keep the reasoning: Modules 3, 4, 5, 6 in order. They are about
transformation bias, resampling, weighted least squares and model selection, using enzyme
kinetics only as a source of realistic data.

### The software route

Skip the statistics, build the tool: read the source of `src/kinetics/`, then
[practicals 1, 4 and 5](practicals/practical-01-first-michaelis-menten-fit.md), then the
capstone. The tests are the specification — start with `tests/test_linearisation.py`,
which encodes every claim the documentation makes about bias.

## Stop points

* **After Module 2** you can fit a clean substrate series and say what you did. That is a
  defensible place to stop if your only need is a Km.
* **After Module 4** you can report uncertainty you would defend in a paper's methods
  section.
* **After Module 6** you can design an experiment, or explain why the one you have cannot
  answer the question. The second is more useful and much rarer.

## What to do when you get stuck

1. Run the command in the module; if it fails, compare your output with the module's. The
   modules print what this repository's code printed on this machine, so a mismatch is
   diagnostic.
2. Read `docs/limitations.md`. Several "errors" are the package refusing to do something
   it cannot do honestly.
3. Run the tests for the area: `python -m pytest tests/test_fitting.py -q`. If they pass,
   the library is fine and the problem is in the data or the setup.
4. If you believe the documentation is wrong, it probably is: open an issue with the
   command and the output. Documentation errors in this repository are treated as defects,
   because a wrong number teaches a wrong habit.
