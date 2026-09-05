---
name: test-suite-shape
description: "How the ported pytest suite is organised, why some files are one big test, and what was deliberately dropped."
metadata:
  node_type: memory
  type: project
  modified: 2026-09-03
---

`python -m pytest` from the repo root. 57 tests, 523 named checks, via a `check(name,
cond, detail)` helper in `code/tests/support.py` -- kept because these assertions describe
a property in words and the name is what tells you which property broke.

## Some files are one test on purpose

The fork's tests were **flat scripts** with `print("[n] label")` section markers, and many
sections share state built by earlier ones. Where sections are independent they were split
into separate test functions; where they are not, the file is **one test function**, which
is the granularity the fork already had. Splitting those would have meant inventing setup
the original never had.

## The two shapes that must both be clean

A checkout **and** a copy with no `.git`. Three provenance defects in earlier runs passed
their tests and failed in the capsule because the tests ran in a shape production never
has. `test_entrypoint` covers the git-less case directly, and the single check that cannot
apply there is `skipif`-marked. Verify with:

```bash
cp -r code docs pytest.ini <tmp>/ && cd <tmp> && python -m pytest -q
```

## Dropped, with reasons

* `test_reference_tables` -- comparison against the historical tables is retired by decision.
* `test_diff_runs`, `test_preflight`, and `test_v1dd_nwb` sections 5-7 -- their targets
  (`compare.py`, `preflight.py`, `schema_report.py`, `checkpoints.py`) were not carried
  into this repo. **If any of those modules is ever ported back, port its test with it.**

`test_entrypoint`, `test_import_boundary` and `test_tuning_export` were rewritten rather
than ported, because the version ladder, the layering and the array writers all changed
shape. A line-by-line port would have tested the old design.

## The ported tests caught five real things

All five were our own deliberate changes, which is what a port is for: the divergence set
growing from five settings to eight, `_ratio` returning NaN instead of 0,
`_lifetime_sparseness_chunked` defaulting to condition means, and the
`surround_supression_index` rename in two places.

One was a **better** catch than that. `test_drifting_gratings` asserted "no dF/F columns"
with the reason *"a ratio index is not safe on a signed trace"* -- a proxy for the real
rule. `run_corr_dff` is a correlation, not a ratio, so it is legitimately on dF/F. The test
now states the actual rule: no *ratio* column on a signed trace, and the one dF/F column
must be the correlation. Fix the rule, not the symptom.
