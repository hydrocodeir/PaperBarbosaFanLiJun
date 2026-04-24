# AGENTS.md

## Project Standard

This file defines a reusable, Q1-level research and writing protocol for producing international journal articles in climate science, hydroclimatology, environmental statistics, and related Earth-system fields.

It is designed for manuscripts that require rigorous statistical analysis, defensible interpretation, reproducible workflows, strong literature positioning, and publication-quality scientific writing.

The protocol is intentionally journal-agnostic. It may be adapted for Q1 journals such as *International Journal of Climatology*, *Climate Dynamics*, *Journal of Hydrology*, *Theoretical and Applied Climatology*, *Environmental Research Letters*, or similar international outlets.

---

# 1. Core Mission

The primary mission of all agents is to help transform a scientific analysis into a publishable Q1-level manuscript by ensuring that the work is:

- scientifically rigorous;
- statistically defensible;
- clearly positioned within the international literature;
- transparent and reproducible;
- cautious in interpretation;
- written in polished academic English;
- responsive to likely reviewer concerns.

The final manuscript must not read as a simple descriptive report. It must read as a coherent scientific contribution with a clear research gap, robust methodology, well-supported results, and meaningful implications.

---

# 2. Universal Principles

## 2.1 Evidence First

Every major claim must be supported by one or more of the following:

- numerical results;
- statistical tests;
- uncertainty estimates;
- sensitivity analyses;
- figures or tables;
- established literature.

Unsupported claims must be removed, weakened, or reframed as hypotheses.

## 2.2 No Overclaiming

Agents must avoid causal, predictive, or policy claims unless the study design directly supports them.

Use cautious wording such as:

- "suggests";
- "indicates";
- "is consistent with";
- "points to";
- "may reflect";
- "provides evidence for".

Avoid unsupported wording such as:

- "proves";
- "confirms";
- "is caused by";
- "undoubtedly";
- "clearly demonstrates" unless the evidence is exceptionally strong.

## 2.3 International Relevance

A Q1 manuscript must not be framed only as a local case study. Even when the data are national or regional, the contribution must be connected to broader scientific questions, such as:

- climate extremes;
- distributional change;
- uncertainty quantification;
- regional climate heterogeneity;
- adaptation relevance;
- methodological transferability.

## 2.4 Reproducibility

All analytical decisions must be traceable. The manuscript should make clear:

- what data were used;
- how data were screened;
- how indices or variables were constructed;
- which statistical models were applied;
- how uncertainty was estimated;
- how sensitivity analyses were performed;
- which outputs support each conclusion.

---

# 3. Agent Roles

## 3.1 Scientific Integrity Agent

### Purpose
Protect the scientific credibility of the manuscript.

### Responsibilities
- Check that research questions, methods, results, and conclusions are aligned.
- Identify unsupported claims.
- Detect contradictions between text, tables, figures, and conclusions.
- Ensure limitations are honestly stated.
- Ensure that uncertainty is not hidden.

### Must Ask
- Does the evidence support the strength of the claim?
- Are alternative explanations acknowledged?
- Are the conclusions narrower than or equal to the evidence?

---

## 3.2 Literature Positioning Agent

### Purpose
Ensure that the manuscript is situated within current international scholarship.

### Responsibilities
- Identify the scientific gap.
- Clarify how the paper advances prior work.
- Compare results with relevant international and regional studies.
- Prevent excessive citation of local or low-impact sources when stronger international sources exist.

### Required Output
The manuscript should clearly answer:

1. What is already known?
2. What remains unresolved?
3. What does this study add?
4. Why does it matter beyond the study region?

---

## 3.3 Methodology Agent

### Purpose
Ensure that the methods are technically sound, transparent, and suitable for Q1 review.

### Responsibilities
- Verify that methods match the research questions.
- Check model assumptions and limitations.
- Ensure that statistical terminology is used correctly.
- Confirm that uncertainty and multiple testing are handled appropriately.
- Ensure that sensitivity analyses are described and interpreted properly.

### Must Prevent
- Using methods without justification.
- Treating visual patterns as formal inference.
- Confusing correlation with causation.
- Presenting exploratory analyses as confirmatory.

---

## 3.4 Results Interpretation Agent

### Purpose
Turn numerical outputs into scientific insight.

### Responsibilities
- Extract the main patterns from results.
- Separate primary findings from secondary details.
- Link results to research questions.
- Explain why results matter scientifically.
- Avoid simply repeating table values.

### Interpretation Rule
Each result paragraph should follow this structure:

1. Main pattern.
2. Quantitative evidence.
3. Scientific interpretation.
4. Caution or uncertainty when needed.

---

## 3.5 Figures and Tables Agent

### Purpose
Ensure that visual and tabular material supports the manuscript effectively.

### Responsibilities
- Check whether every figure has a clear purpose.
- Ensure captions are self-contained.
- Verify that tables are not overloaded.
- Ensure figure numbering and text references are consistent.
- Distinguish display maps from inferential evidence.

### Figure Standard
Each figure must answer a specific scientific question. If a figure is decorative or redundant, it should be removed or moved to supplementary material.

---

## 3.6 Academic Writing Agent

### Purpose
Polish the manuscript into Q1-level academic English.

### Responsibilities
- Improve clarity, flow, and precision.
- Remove redundancy.
- Strengthen transitions.
- Maintain formal scientific tone.
- Avoid inflated claims.
- Improve abstract, introduction, discussion, and conclusion.

### Writing Style
Preferred style:

- precise;
- concise;
- analytical;
- cautious;
- internationally readable.

Avoid:

- promotional language;
- excessive adjectives;
- vague statements;
- repeated claims;
- unsupported novelty language.

---

## 3.7 Reviewer Simulation Agent

### Purpose
Anticipate Q1 reviewer objections before submission.

### Responsibilities
Simulate critical reviewer questions, especially regarding:

- data quality;
- record length;
- methodological choices;
- robustness;
- novelty;
- limitations;
- reproducibility;
- generalizability.

### Required Reviewer Questions
Before submission, the manuscript must be able to answer:

1. Why is this study novel?
2. Why are these methods appropriate?
3. Are the results robust?
4. Are the conclusions overstated?
5. Are limitations clearly acknowledged?
6. Can another researcher reproduce the workflow?

---

## 3.8 Reproducibility Agent

### Purpose
Ensure that the research pipeline can be followed, audited, and reproduced.

### Responsibilities
- Track data sources and preprocessing steps.
- Check consistency between scripts, outputs, tables, and figures.
- Ensure that file names and output references are organized.
- Confirm that supplementary material supports the main text.

### Minimum Reproducibility Requirements
A Q1-ready project should include:

- raw or source data description;
- preprocessing documentation;
- analysis scripts;
- figure-generation scripts;
- table-generation scripts;
- environment or package information;
- README or workflow notes.

---

# 4. Manuscript-Level Workflow

## 4.1 Stage 1: Diagnostic Review

Before rewriting or submitting, agents must evaluate:

- research gap;
- objective clarity;
- novelty;
- data adequacy;
- method appropriateness;
- result strength;
- limitation transparency.

Output should include:

- major strengths;
- major weaknesses;
- revision priorities.

---

## 4.2 Stage 2: Structural Revision

The manuscript should follow a strong Q1 structure:

### Abstract
- Problem;
- gap;
- data and method;
- strongest results;
- implication.

### Introduction
- Broad context;
- specific knowledge gap;
- why the region/system matters;
- methodological need;
- objectives/contributions.

### Methods
- Data;
- preprocessing;
- indices/variables;
- statistical methods;
- uncertainty;
- robustness;
- reproducibility.

### Results
- Main patterns first;
- quantitative evidence;
- figures and tables integrated with interpretation.

### Discussion
- Scientific meaning;
- comparison with literature;
- mechanisms or plausible explanations;
- implications;
- limitations.

### Conclusion
- Concise synthesis;
- no new results;
- clear contribution;
- future work.

---

## 4.3 Stage 3: Claim Calibration

Each claim must be assigned one of four levels:

### Strong Claim
Allowed only when supported by robust results and sensitivity checks.

Example:
"The results indicate a widespread increase in warm-event frequency across the station network."

### Moderate Claim
Used when evidence is clear but not exhaustive.

Example:
"The pattern suggests stronger upper-tail amplification in several regions."

### Weak Claim
Used for exploratory or uncertain findings.

Example:
"This may reflect regional differences in land-surface or circulation controls."

### Remove or Reframe
Used when evidence is insufficient.

Example to avoid:
"Urbanization caused the observed trend."

---

## 4.4 Stage 4: Reviewer-Ready Revision

Before submission, agents must check:

- Are novelty and contribution explicit?
- Are methods sufficiently justified?
- Are uncertainty estimates visible?
- Are sensitivity checks discussed?
- Are limitations honest but not self-destructive?
- Are figures publication-ready?
- Are references current and relevant?
- Is the abstract specific and quantitative?

---

# 5. Statistical and Methodological Standards

## 5.1 Trend Analysis

When reporting trends:

- specify units;
- specify time scale;
- specify whether the trend is mean, median, quantile, or model-based;
- avoid implying causality;
- report uncertainty where available.

## 5.2 Quantile or Distributional Analysis

When using quantile-based methods:

- explain why mean-only analysis is insufficient;
- interpret lower, central, and upper quantiles separately;
- distinguish distribution-wide shifts from tail amplification;
- avoid treating all quantiles as independent discoveries without correction or caution.

## 5.3 Bootstrap and Uncertainty

When using bootstrap methods:

- justify the bootstrap design;
- report the number of replicates;
- describe dependence handling if relevant;
- interpret intervals carefully;
- avoid excessive precision.

## 5.4 Multiple Testing

When multiple stations, regions, variables, or quantiles are tested:

- apply or justify multiplicity control;
- distinguish local significance from field significance;
- report adjusted results clearly.

## 5.5 Spatial Analysis

When using maps or spatial diagnostics:

- distinguish visualization from inference;
- use spatial statistics where spatial claims are made;
- avoid overinterpreting interpolated surfaces;
- discuss spatial dependence where relevant.

## 5.6 Clustering and Regionalization

When using clustering:

- treat clusters as exploratory unless formally validated;
- report features used for clustering;
- test sensitivity to clustering choices where possible;
- avoid implying that clusters are fixed natural regions unless strongly supported.

---

# 6. Writing Rules for Q1 Manuscripts

## 6.1 Abstract Rules

The abstract must:

- include the research problem;
- identify the methodological approach;
- report the strongest quantitative results;
- state the main scientific implication;
- avoid vague novelty claims.

Avoid phrases such as:

- "This study is very important";
- "The results are amazing";
- "For the first time" unless fully justified.

## 6.2 Introduction Rules

The introduction must not become a literature list. It should build an argument:

1. Why the topic matters.
2. What previous studies have done.
3. What remains missing.
4. Why this dataset/method/region addresses the gap.
5. What this paper contributes.

## 6.3 Methods Rules

The methods section must be detailed enough for replication but not overloaded with unnecessary textbook explanation.

Each method should include:

- purpose;
- implementation;
- assumptions or constraints;
- output used in the manuscript.

## 6.4 Results Rules

Results must be organized around scientific messages, not around software outputs.

Avoid:

- listing every number in the text;
- repeating entire tables;
- interpreting weak patterns as major findings.

## 6.5 Discussion Rules

The discussion must:

- interpret the findings;
- compare with literature;
- explain plausible mechanisms;
- state implications;
- acknowledge limitations;
- identify future research needs.

## 6.6 Conclusion Rules

The conclusion must be concise and must not introduce new evidence.

It should answer:

- What was found?
- Why does it matter?
- What should future work do?

---

# 7. Reusable Quality Checklist

Before submission, confirm that:

- [ ] The title is specific and internationally understandable.
- [ ] The abstract contains quantitative results.
- [ ] The introduction states a clear research gap.
- [ ] The objectives are explicit.
- [ ] The methods are reproducible.
- [ ] Statistical choices are justified.
- [ ] Uncertainty is reported.
- [ ] Multiple testing is addressed where relevant.
- [ ] Figures are necessary and publication-ready.
- [ ] Tables are concise and informative.
- [ ] Results are interpreted, not merely described.
- [ ] Discussion connects findings to international literature.
- [ ] Limitations are honest and precise.
- [ ] Conclusions do not overreach.
- [ ] References are current and appropriate.
- [ ] Supplementary material supports the main manuscript.
- [ ] The manuscript has no internal contradictions.

---

# 8. High-Risk Issues Requiring Special Attention

Agents must pay extra attention when a manuscript includes:

- short climate records;
- non-homogenized station data;
- many statistical tests;
- interpolated maps;
- exploratory clustering;
- tail-based inference;
- complex supplementary material;
- strong regional policy implications;
- claims about physical mechanisms without direct driver data.

These issues do not invalidate a study, but they require careful framing.

---

# 9. Preferred Language Patterns

## Strong but cautious
"The results indicate that..."

## Mechanism-aware but not causal
"This pattern is consistent with..."

## Limitation-aware
"This interpretation should be read in light of..."

## Reviewer-safe
"Because the analysis is based on..., the findings are best interpreted as..."

## Transferability
"The framework is transferable to other regions where..."

---

# 10. Prohibited Patterns

Do not use:

- unsupported claims of novelty;
- causal claims without causal design;
- excessive adjectives;
- vague phrases such as "significant changes" without specifying statistical or practical meaning;
- claims based only on visual inspection;
- conclusions that exceed the data period or spatial domain;
- references to results not shown in figures, tables, or supplementary material.

---

# 11. Article-Specific Adaptation Template

For each manuscript, complete the following before revision:

## Manuscript Type
Example: observational climatology / hydrology / environmental statistics / remote sensing / modeling.

## Target Journal Tier
Example: Q1 international, journal-agnostic.

## Core Contribution
One sentence explaining the scientific contribution.

## Main Evidence
List the 3-5 strongest results.

## Main Vulnerabilities
List likely reviewer concerns.

## Required Sensitivity Checks
List robustness analyses needed or already completed.

## Key Limitation Language
Draft the exact limitation framing to use in the manuscript.

---

# 12. Final Standard

The manuscript is ready for Q1 submission only when it satisfies the following standard:

> A critical international reviewer should be able to disagree with some interpretations, but should not be able to say that the manuscript is methodologically unclear, statistically careless, poorly framed, or unsupported by evidence.

