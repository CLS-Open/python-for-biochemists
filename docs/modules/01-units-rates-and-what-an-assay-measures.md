# Module 1 — Units, rates and what an assay actually measures

**Time:** 75 minutes. **Prerequisite:** none beyond [installation](../installation.md).

## Learning objectives

By the end of this module you can:

1. state an initial rate completely — numerator, denominator, temperature and buffer — and say
   why a rate without a denominator is not a result;
2. convert between µmol min⁻¹ mg⁻¹, mol s⁻¹ mg⁻¹ and the same quantity per unit of enzyme;
3. explain why two Km values with different denominators can be compared but two Vmax values
   cannot, and what else has to match for either comparison to mean anything;
4. read a kinetics CSV before fitting it, and name the four things you check first (units,
   concentration range, replicate structure, obviously impossible values).

---

## 1.1 Why this module exists

Ask ten people for "the Km of hexokinase" and you will get ten numbers with no units, from
0.02 to 10 mM, all of them defensible. The differences come from what was measured, how
the rate was defined, what the denominator was, and at what pH and temperature. None of
those is a property of the enzyme alone.

A kinetics number without its measurement conditions is not a measurement. It is a
recollection. This module is about the conditions, and the rest of the course is only
meaningful if you can state them.

## 1.2 The rate is a slope, not a point measurement

An initial rate is `d[P]/dt` in the early part of the reaction, where the product-versus-time
curve is still straight. The reason is not convenience:

* as substrate is consumed, `[S]` falls, so the rate falls even if the enzyme is unchanged;
* product often inhibits, so the rate falls faster than substrate depletion alone predicts;
* the enzyme may lose activity at 37 °C over tens of minutes.

The straight part of the curve is the only part where the measured d[P]/dt can be assigned
to the `[S]` you pipetted. Everything after it is a mixture of effects.

> **Decision you are making when you choose the linear window:** you are asserting that
> over that window, substrate depletion and product inhibition are negligible. Write down
> the window. A rate measured over 5 minutes from a reaction with Km = 20 µM and
> S₀ = 20 µM has consumed about 30 % of the substrate and is not an initial rate.

## 1.3 The denominator problem

A rate in `µmol min⁻¹` is not comparable between two experiments, because it depends on how
much enzyme was in the cuvette. Kinetic constants are reported per amount of enzyme:

| Denominator | Typical unit | When it is right |
|---|---|---|
| Total protein | µmol min⁻¹ (mg protein)⁻¹ | crude lysate, where you cannot see the enzyme |
| Purified enzyme mass | µmol min⁻¹ mg⁻¹ | pure enzyme, specific activity |
| Enzyme active sites | µmol min⁻¹ (µmol sites)⁻¹ = min⁻¹ | `kcat`, the only denominator that is a molecular property |
| Cell number | µmol min⁻¹ (10⁶ cells)⁻¹ | whole-cell assays |
| Reaction volume | µmol min⁻¹ mL⁻¹ | almost never useful; it measures your cuvette |

`Vmax` inherits the denominator. `Km` does **not** — it is a concentration, so it is
independent of how much enzyme you used, which is exactly why it is the more portable
number and exactly why it is tempting to over-interpret.

**This course uses `µmol min⁻¹ (mg protein)⁻¹` for rates and `µM` for concentrations**,
and the package never converts units. If you feed it mM, you get Km in mM and the same
Vmax in your rate units. That is not laziness: automatic conversion requires guessing
whether `25` means µM or mM, and a guess that is wrong by 1000× is worse than no answer.
[Module 2](02-the-michaelis-menten-model-and-the-two-numbers-people-misquote.md) shows a
test that pins this behaviour.

## 1.4 Temperature, pH, buffer, ionic strength

Every one of these changes Vmax and often Km. The classic example is that Km for a
substrate can shift several-fold over a pH unit if the binding group titrates in that
range. So:

* state temperature and pH next to the constant;
* do not compare constants measured at different conditions and call the difference a
  mechanism;
* if you must compare across papers, compare *within* one set of conditions and say so.

A Km quoted at "25 °C, pH 7.4, 50 mM Tris-HCl, 5 mM MgCl₂" is a measurement. The same
number quoted alone is a rumour with a decimal point.

## 1.5 What the substrate concentration actually is

`[S]` in the rate law is the **free** substrate concentration at steady state, not the
amount you pipetted. They differ when:

* the substrate is substantially consumed by the assay (see 1.2);
* the substrate binds the enzyme at a stoichiometric ratio (tight-binding inhibition, and
  `[E]` comparable to `Ki` — outside the scope of this course, see
  [limitations](../limitations.md));
* the substrate is buffered by a chelator: `Mg·ATP` is the real substrate for kinases, not
  `ATP`.

The last one is the most common silent error in student data: an ATP titration done in
excess EDTA, or with `Mg²⁺` held below the ATP concentration, produces a Km that is
reproducible in that buffer and meaningless outside it.

## 1.6 Looking at data before fitting it

Rule: plot the raw data before any model touches it. The repository gives you a way to
check that a file is even the right shape:

```bash
kinetics mm datasets/sample/malformed_example.csv; echo "exit=$?"
```

```
error: datasets/sample/malformed_example.csv line 3: v='3.4x' is not a number
exit=2
```

Exit code 2 means "this file cannot be read", and the file is in the repository to be
rejected: it contains a non-numeric cell and a negative rate. A tool that quietly dropped
those rows would give you a fit and no warning.

The other datasets are read this way — no library, no `pandas`, just the format:

```python
import csv
from pathlib import Path

rows = list(csv.DictReader(Path("datasets/sample/michaelis_menten.csv").open(encoding="utf-8")))
print(len(rows), "rows;", list(rows[0]))
print(rows[:3])
```

```
27 rows; ['S', 'replicate', 'v']
[{'S': '1.0', 'replicate': '1', 'v': '0.33515'},
 {'S': '1.0', 'replicate': '2', 'v': '0.37198'},
 {'S': '1.0', 'replicate': '3', 'v': '0.37316'}]
```

Three replicates at each of nine substrate concentrations. Every value is a string at this
point; the conversion to `float` is a decision, and the CLI refuses rather than guesses
when it cannot do it.

## 1.7 The units table you should copy into your notebook

```python
UNIT_TABLE = {
    "S": "µM",  # substrate, free concentration
    "I": "µM",  # inhibitor
    "v": "µmol min^-1 (mg protein)^-1",  # initial rate
    "vmax": "same as v",
    "km": "µM",
    "ki": "µM",
    "kcat": "min^-1",  # per active site, needs [E]
}
```

Two lines in a notebook cost nothing and prevent the single most common review comment on
kinetics manuscripts.

## 1.8 Worked check

`datasets/sample/michaelis_menten.csv` was generated with `Vmax = 10.0` and `Km = 25.0`
µM. At `[S] = 20 µM` the true rate is therefore `10·20/45 = 4.44`. The measured values at
that concentration are 4.28, 4.39 and 4.37:

```python
import csv

rows = list(csv.DictReader(open("datasets/sample/michaelis_menten.csv", encoding="utf-8")))
print([r["v"] for r in rows if r["S"] == "20.0"])
```

```
['4.28471', '4.39016', '4.37152']
```

If the three values disagree by more than you would accept in your own assay, the
simulated noise level is wrong for your purposes — and that is a *useful* complaint,
because you now know what the file contains rather than assuming.

## 1.9 Exercises

1. The three rates measured at `[S] = 100 µM` are 7.63, 8.25 and 8.02. Compute the true
   value from `Vmax` and `Km` and check that the mean is within the noise. What fraction of
   Vmax is reached at 100 µM, and what would you need to measure to see 95 % of it?
2. Suppose you mistakenly divided all rates by 2 (forgot half the enzyme). Which of Vmax,
   Km, and the fitted R² change?
3. A colleague reports "Km = 0.02" with no units. List the three questions you need to ask
   before you can use it.
4. Why does a rate measured at `[S] ≈ Km` over 10 % substrate depletion bias Km downward
   more than upward?

## 1.10 What this module did **not** cover

* Steady-state derivation of the rate law (Module 2 sketches it; a full derivation needs
  the Briggs–Haldane treatment).
* Tight-binding inhibition, where `[E]` is not negligible against `[I]` (see
  [limitations](../limitations.md)).
* Coupled assays, where the rate you measure is another enzyme's.
* Stopped-flow and pre-steady-state kinetics: this repository is about steady-state
  initial rates only.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `kinetics: error: unrecognized arguments: --weighted` | There is no such flag: weighting is triggered by the presence of an `sd` column | Add the column, or fit unweighted and say so |
| `error: datasets/sample/... : no such file` | You are not in the repository root | `cd` to the repository root; every documented path is relative to it |
| A rate that looks 1000× too small | µmol vs mol, or per mL vs per mg | Write the conversion out; the tool does no unit conversion and will not warn you |
| Two columns of numbers, no units anywhere | Units live in the metadata, not the CSV | Read `datasets/metadata/*.yml`; put the units in your methods sentence |

---

**Next:** [Module 2 — the Michaelis–Menten model and the two numbers people misquote](02-the-michaelis-menten-model-and-the-two-numbers-people-misquote.md)
