## What does this PR change?

<!-- One paragraph. Link the issue it closes: "Closes #123" -->

## Type

- [ ] Content correction (scientific, factual, typo)
- [ ] New or revised lesson / practical / assessment
- [ ] Software change
- [ ] Documentation / infrastructure
- [ ] Translation

## Scientific-accuracy checklist (required for content PRs)

- [ ] Every new factual claim has a reference that I opened and read
- [ ] Every new reference resolves via Crossref/Europe PMC (checked on: ________)
- [ ] Arithmetic in examples verified independently (how? ________)
- [ ] Units stated for every number
- [ ] Synthetic/simulated data is labelled as such in text **and** in figures
- [ ] Limitations of the method taught are stated
- [ ] No fabricated data, citations, authors or results

## Software checklist (required for code PRs)

- [ ] `pytest -q` passes locally
- [ ] New logic has tests; arithmetic expectations are hand-computed in comments
- [ ] Public functions have docstrings stating complexity/limitations where relevant
- [ ] No new dependency unless justified in the PR description
- [ ] No secrets, tokens or personal data added anywhere

## Educational checklist

- [ ] Learning objectives are measurable and stated at the top
- [ ] Exercise has a deliverable and a rubric
- [ ] Difficulty matches the module level
- [ ] Accessibility: colour-blind-safe colours, alt text, no reliance on colour alone

## Evidence of reproducibility

<!-- Paste the commands you ran and their result, e.g. `python -m pytest -q` → 230 passed -->

## Anything you are unsure about?

<!-- Reviewers would rather hear it from you than find it later. -->
