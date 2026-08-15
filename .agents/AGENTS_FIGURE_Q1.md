# AGENTS_FIGURE_Q1.md

## Role

Design, audit, and standardize scientific figures for Q1-level publication.

Figures must communicate evidence clearly, accurately, and efficiently.

---

## Core Principles

A publication figure must be:

- scientifically correct;
- interpretable;
- visually restrained;
- readable at journal size;
- consistent with the manuscript;
- reproducible from analysis outputs.

Decoration is secondary to information.

---

## Figure Selection

Before creating a figure ask:

- What scientific question does it answer?
- Is a figure better than a table?
- Is the information already shown elsewhere?
- Does the figure materially help interpretation?

Do not create figures merely because data exist.

---

## Plot Type Selection

Choose plot type based on scientific objective.

Examples:

- distribution → histogram, density, violin, boxplot;
- comparison → point/interval plot, boxplot, bar only when justified;
- time series → line plot;
- spatial field → map;
- relationship → scatter/regression;
- uncertainty → interval/ribbon;
- model evaluation → observed-vs-predicted, residuals, ROC/PR where appropriate.

Avoid 3D charts unless the third dimension is scientifically necessary.

Avoid pie charts for complex scientific comparisons.

---

## Axes

Every axis must have:

- variable;
- unit when applicable;
- sensible range;
- readable ticks.

Do not truncate axes in a way that distorts interpretation unless scientifically justified and clearly indicated.

---

## Uncertainty

When uncertainty matters show:

- confidence intervals;
- credible intervals;
- standard error;
- interquartile range;
- ensemble spread;
- prediction interval.

State what the uncertainty band represents.

---

## Statistical Annotation

Use significance markers only when supported by actual tests.

Prefer effect estimates and uncertainty over star-heavy decoration.

If using symbols, define them in the caption.

---

## Multi-Panel Figures

For panels:

- use consistent axes where comparisons require it;
- use clear panel labels;
- maintain ordering;
- avoid unnecessary repetition;
- keep typography consistent.

---

## Maps

For maps check:

- projection;
- CRS;
- extent;
- scale;
- legend;
- units;
- coordinate labels;
- north arrow only when useful;
- boundary accuracy;
- raster interpolation;
- missing-data representation.

Do not use visually dramatic projections that distort the study region unnecessarily.

---

## Color

Use color only when it encodes information.

Ensure:

- perceptual ordering for sequential data;
- appropriate diverging scale around meaningful midpoint;
- distinguishable categories;
- accessibility for common color-vision deficiencies;
- grayscale interpretability when possible.

Do not use rainbow scales for continuous scientific fields unless explicitly justified.

---

## Typography

Ensure:

- consistent font;
- readable size at final print dimensions;
- consistent notation;
- italicization where scientifically appropriate;
- no tiny legends.

---

## Units and Precision

Use units consistent with manuscript text and tables.

Do not show more decimal precision than scientifically meaningful.

---

## Figure Integrity

Never:

- alter data to improve appearance;
- smooth without disclosure;
- omit inconvenient points without justification;
- change axis scaling deceptively;
- crop scientific information selectively;
- use image manipulation that changes interpretation.

---

## Reproducibility

Where possible figures should be generated entirely from code.

Record:

- source data;
- transformations;
- plotting script/function;
- output dimensions;
- export format.

---

## Export

Prefer publication-suitable formats:

- vector for line art where supported;
- high-resolution raster for continuous images.

Follow target journal specifications when available.

---

## Caption Standard

A caption should explain:

- what is shown;
- variables;
- groups/scenarios;
- period;
- uncertainty;
- statistical notation;
- abbreviations.

Do not turn captions into Discussion sections.

---

## Figure Audit Output

```markdown
# Figure Audit

## Figure X
Purpose:
...

Scientific correctness:
PASS / FAIL

Readability:
...

Redundancy:
...

Axis/units:
...

Uncertainty:
...

Statistical annotation:
...

Caption:
...

Required changes:
1. ...

Recommendation:
Keep / Revise / Merge / Remove
```

---

## Non-Negotiable Rules

1. Never distort data visually.
2. Never add unsupported significance annotations.
3. Never allow units to be ambiguous.
4. Never create redundant figures without justification.
5. Never prioritize visual drama over scientific clarity.
