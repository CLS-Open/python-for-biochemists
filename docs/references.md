# References

Every DOI below was resolved by `scripts/check_references.py`, which compares the title
printed here with the title Crossref returns for that DOI. A DOI that resolves to a different
paper is a worse failure than a dead link, because it looks correct; the script exists because
that happened during the writing of this repository.

```bash
python scripts/check_references.py            # resolves every DOI via Crossref
python scripts/check_references.py --offline  # internal consistency only
```

## The model and its two parameters

| Work | Title | DOI | Why it is here |
|---|---|---|---|
| Hill (1913) | **The Combinations of Haemoglobin with Oxygen and with Carbon Monoxide. I** | [10.1042/bj0070471](https://doi.org/10.1042/bj0070471) | the Hill equation, and where the shape parameter comes from |
| Hanes (1932) | **Studies on plant amylases** | [10.1042/bj0261406](https://doi.org/10.1042/bj0261406) | Hanes–Woolf: the first linearisation, and the one that damages least |
| Lineweaver & Burk (1934) | **The Determination of Enzyme Dissociation Constants** | [10.1021/ja01318a036](https://doi.org/10.1021/ja01318a036) | Lineweaver–Burk: the double reciprocal this repository measures rather than repeats |
| HOFSTEE (1959) | **Non-Inverted Versus Inverted Plots in Enzyme Kinetics** | [10.1038/1841296b0](https://doi.org/10.1038/1841296b0) | Hofstee's rebuttal, and the basis of the Eadie–Hofstee diagnostic |
| Cleland (1963) | **The kinetics of enzyme-catalyzed reactions with two or more substrates or products** | <https://doi.org/10.1016/0006-3002(63)91800-6> | Cleland's nomenclature: competitive, uncompetitive, noncompetitive |

## Fitting, weighting and uncertainty

| Work | Title | DOI | Why it is here |
|---|---|---|---|
| Marquardt (1963) | **An Algorithm for Least-Squares Estimation of Nonlinear Parameters** | [10.1137/0111030](https://doi.org/10.1137/0111030) | the Levenberg–Marquardt algorithm that SciPy's curve_fit uses |
| Nelder & Mead (1965) | **A Simplex Method for Function Minimization** | [10.1093/comjnl/7.4.308](https://doi.org/10.1093/comjnl/7.4.308) | Nelder–Mead, the optimiser you fall back to when derivatives are unavailable |
| WILKINSON (1961) | **Statistical estimations in enzyme kinetics** | [10.1042/bj0800324](https://doi.org/10.1042/bj0800324) | Wilkinson: why reciprocal plots need weights, written when it still had to be argued |
| Cornish-Bowden & Eisenthal (1974) | **Statistical considerations in the estimation of enzyme kinetic parameters by the direct linear plot and other methods** | [10.1042/bj1390721](https://doi.org/10.1042/bj1390721) | Cornish-Bowden and Eisenthal: direct linear plots, and the statistics of the alternatives |
| Rocke & Lorenzato (1995) | **A Two-Component Model for Measurement Error in Analytical Chemistry** | [10.1080/00401706.1995.10484302](https://doi.org/10.1080/00401706.1995.10484302) | the two-component error model used by every dataset in this repository |
| Hurlbert (1984) | **Pseudoreplication and the Design of Ecological Field Experiments** | [10.2307/1942661](https://doi.org/10.2307/1942661) | Hurlbert on pseudoreplication: why three replicates are not 27 points |
| Martin Bland & Altman (1986) | **STATISTICAL METHODS FOR ASSESSING AGREEMENT BETWEEN TWO METHODS OF CLINICAL MEASUREMENT** | <https://doi.org/10.1016/S0140-6736(86)90837-8> | Bland–Altman: how to compare two methods without a correlation coefficient |

## Model choice and resampling

| Work | Title | DOI | Why it is here |
|---|---|---|---|
| Akaike (1974) | **A new look at the statistical model identification** | [10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705) | Akaike's criterion |
| HURVICH & TSAI (1989) | **Regression and time series model selection in small samples** | [10.1093/biomet/76.2.297](https://doi.org/10.1093/biomet/76.2.297) | Hurvich and Tsai: the small-sample correction without which AIC picks the wrong model here |
| Burnham & Anderson (2004) | **Multimodel Inference** | [10.1177/0049124104268644](https://doi.org/10.1177/0049124104268644) | Burnham and Anderson on multimodel inference and what ΔAICc does and does not mean |
| Efron (1979) | **Bootstrap Methods: Another Look at the Jackknife** | [10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552) | Efron's bootstrap |
| DiCiccio & Efron (1996) | **Bootstrap confidence intervals** | [10.1214/ss/1032280214](https://doi.org/10.1214/ss/1032280214) | DiCiccio and Efron on bootstrap confidence intervals, including where percentile intervals fail |

## Calibration, detection and quantification limits

| Work | Title | DOI | Why it is here |
|---|---|---|---|
| Currie (1968) | **Limits for qualitative detection and quantitative determination. Application to radiochemistry** | [10.1021/ac60259a007](https://doi.org/10.1021/ac60259a007) | Currie's original detection-limit framework |
| Currie (1995) | **Nomenclature in evaluation of analytical methods including detection and quantification capabilities (IUPAC Recommendations 1995)** | [10.1351/pac199567101699](https://doi.org/10.1351/pac199567101699) | the IUPAC nomenclature for LOD and LOQ |
| Nash (1953) | **The colorimetric estimation of formaldehyde by means of the Hantzsch reaction** | [10.1042/bj0550416](https://doi.org/10.1042/bj0550416) | Nash: the colorimetric chemistry behind the calibration dataset's shape |

## Reproducibility and software

| Work | Title | DOI | Why it is here |
|---|---|---|---|
| Sandve et al. (2013) | **Ten Simple Rules for Reproducible Computational Research** | [10.1371/journal.pcbi.1003285](https://doi.org/10.1371/journal.pcbi.1003285) | ten rules this repository tries to follow, including 'everything from code' |
| Virtanen et al. (2020) | **SciPy 1.0: fundamental algorithms for scientific computing in Python** | [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) | SciPy, whose curve_fit and stats functions every number here is checked against |

## Primary sources without a DOI

Two papers that everything here descends from were published before DOIs existed and have no
registrar record that resolves cleanly. They are listed without links rather than with a DOI
that belongs to something else:

* **Michaelis, L. & Menten, M. L. (1913).** Die Kinetik der Invertinwirkung. *Biochemische
  Zeitschrift* **49**, 333–369. The original derivation, and the source of the Vmax/Km
  parameterisation this repository uses.
* **Hill, A. V. (1910).** The possible effects of the aggregation of the molecules of
  haemoglobin on its dissociation curves. *Journal of Physiology* **40** (suppl.), iv–vii.
  The proceedings abstract that introduced the Hill equation. The full paper with a DOI is
  the 1913 *Biochemical Journal* article listed above.

## Textbooks, not cited by DOI

Books are not checkable by DOI and are therefore not in the table. The two that shaped the
code are:

* Cornish-Bowden, A. *Fundamentals of Enzyme Kinetics*, 4th ed. — the reference for the
  inhibition nomenclature and for why the direct linear plot exists.
* Motulsky, H. & Christopoulos, A. *Fitting Models to Biological Data Using Linear and
  Nonlinear Regression* — the practical source for the AICc and weighting advice in Module 6.

