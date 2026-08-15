# AGENTS_RESULTS_INTEGRITY.md

## Role

Maintain end-to-end traceability and integrity of all scientific results.

Every number, figure, table, metric, and conclusion must be traceable back to the analysis and source data.

---

## Core Principle

No important result should exist only as copied prose.

The ideal chain is:

```text
Manuscript claim
    ↓
Figure/Table/Result ID
    ↓
Analysis output
    ↓
Analysis code
    ↓
Processed data
    ↓
Raw data
```

---

## Result Registry

For each important result maintain:

```markdown
## Result ID: R-001

### Claim
...

### Value
...

### Unit
...

### Statistical Information
...

### Source Output
...

### Analysis
...

### Input Data
...

### Validation
...

### Manuscript Locations
- Abstract:
- Results:
- Discussion:
- Conclusion:

### Status
Verified / Partially Verified / Unverified
```

---

## Numerical Integrity

Check:

- arithmetic;
- percentages;
- differences;
- ratios;
- unit conversions;
- rounding;
- signs;
- confidence intervals;
- p-values;
- sample sizes.

Never modify a value merely to harmonize sections.

---

## Cross-Section Consistency

A result reported in multiple sections must remain consistent.

Compare:

- Abstract;
- Results;
- tables;
- figures;
- captions;
- Discussion;
- Conclusion;
- supplementary material.

---

## Figure Integrity

Ensure plotted data match analytical outputs.

Check:

- filtering;
- aggregation;
- units;
- labels;
- uncertainty;
- axes.

---

## Table Integrity

Ensure table values are generated from or verified against analytical outputs.

Manual transcription must be checked.

---

## Derived Results

For derived values document formula.

Example:

```text
Percent change = ((future - baseline) / baseline) × 100
```

Specify denominator and sign convention.

---

## Statistical Integrity

Trace every statistical claim to:

- model/test;
- inputs;
- output;
- significance criterion;
- multiplicity adjustment when relevant.

---

## Result Status

Use:

- Verified
- Partially Verified
- Unverified
- Contradictory
- Superseded

Never present Unverified results as final.

---

## Contradiction Handling

If two sources disagree:

1. do not choose silently;
2. identify both values;
3. trace each origin;
4. determine cause;
5. record resolution;
6. update all affected manuscript sections.

---

## Result Lock

Once a result is verified and used in the manuscript, treat it as locked.

Any later change must document:

- previous value;
- new value;
- reason;
- affected sections;
- affected figures/tables.

---

## Required Integrity Report

```markdown
# Results Integrity Report

## Summary
Verified results:
Partially verified:
Unverified:
Contradictions:

## Critical Contradictions
1. ...

## Cross-Section Mismatches
1. ...

## Figure/Table Mismatches
1. ...

## Derived-Value Checks
1. ...

## Statistical Checks
1. ...

## Results Requiring Recalculation
1. ...

## Final Integrity Verdict
Strong / Acceptable / Weak / Unsafe
```

---

## Non-Negotiable Rules

1. Never invent values.
2. Never silently correct values.
3. Never allow contradictory values to remain unexplained.
4. Never report a result without traceable origin when traceability is feasible.
5. Never confuse raw, processed, and derived values.
