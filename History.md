# Project history

This file is an append-only record of substantive manuscript, analysis, and reproducibility work.

## 2026-08-15 — Full scientific audit and reference-based rewrite initiated

- User requested a strict-review audit of `reports/Manuscript.md`, all analysis code and outputs, addition of scientifically necessary analyses, and a major Q1-level rewrite saved under a new name.
- User identified `.agents/AGENTS_AI_WRITING.md` as a mandatory writing standard and requested `.agents/AGENTS_Reference-Based-Academic-Rewriter.md` with `assets/Paper 3.pdf` as the style benchmark.
- Existing uncommitted repository changes were detected in the manuscript and analysis code. They are treated as user-owned and will not be reset or overwritten. The revised manuscript will be written to a new file.
- Initial inventory found 124-station daily data for 1991–2024, a configuration-driven Python pipeline, 610 generated outputs, and dedicated thermal-extreme, climate-fingerprint, climate-regime, robustness, and compound dry-hot outputs.
- The reference paper was identified as Alizadeh et al. (2020), *Science Advances*, “A century of observations reveals increasing likelihood of continental-scale compound dry-hot extremes,” DOI 10.1126/sciadv.aaz4571.
- Initial reviewer-risk findings: the current title and abstract are overextended; Results and Discussion are combined and repetitive; several physical interpretations exceed directly tested drivers; the 34-year record requires careful qualification of tail and empirical return-period results; and manuscript wording describes the configured ±5-day threshold window ambiguously as a “5-day moving window” although the implementation uses an 11-day window centered on each calendar day.
- Next steps: complete code/output integrity checks, verify critical citations and numerical claims, decide whether an additional analysis is necessary, implement and rerun only justified additions, then produce and audit a newly named manuscript.

## 2026-08-15 — Additional analyses selected and implemented

- Confirmed 45 daily records with Tmin greater than Tmax (32 stations). A second, non-overlapping screen identified 2,165 records with Tmean outside the valid Tmin–Tmax interval (61 stations). The earlier manuscript count of 2,336/78 included the Tmin–Tmax conflict rows and therefore did not represent mutually exclusive screening categories.
- Added a non-destructive internal-consistency sensitivity analysis. The source observations and primary outputs remain unchanged; the sensitivity masks contradictory Tmin/Tmax/Tmean values, rebuilds the annual indices, and compares regional quantile slopes and compound dry-hot summaries.
- The screening had negligible effects on the regional thermal results: the largest absolute change among the reported network slopes was 0.026 days per decade, and the largest change in Delta1 was 0.024 days per decade. Compound affected-station means changed by at most 0.53 station and linear trend slopes by at most 0.38 station per decade.
- Added a residual moving-block sensitivity for compound affected-station trends (4-year blocks, 4,999 replicates). All reported RP >= 10 and RP >= 20 trends retained positive dependence-aware evidence; for warm-season events, block-null p = 0.0002 for both Kendall tau and linear slope at both thresholds. Annual-event Kendall p-values were 0.0152 (RP >= 10) and 0.0258 (RP >= 20).
- Added output tables: `temperature_internal_consistency_screening.csv`, `temperature_internal_consistency_quantile_sensitivity.csv`, `temperature_internal_consistency_compound_sensitivity.csv`, and `compound_dry_hot_serial_dependence_sensitivity.csv`.
- Added `requirements.txt` after discovering that the repository did not document runtime dependencies. The current environment lacks statsmodels, so the full historical pipeline cannot be rerun in this environment without installing that declared dependency. Targeted new analyses were executed using installed scientific packages while preserving the production pipeline's statsmodels requirement.
- Verified the existence and bibliographic details of several recent references used in the manuscript, including Jamali et al. (2026), Schuhen et al. (2026), Peng et al. (2026), Perkins-Kirkpatrick et al. (2026), and Lv et al. (2026).

## 2026-08-15 — Q1 rewrite, supplementary index, and final audit completed

- Created `reports/Manuscript_Q1_Revised.md` under a new name, preserving the user's original manuscript. The new title is “Quantile-dependent changes in thermal extremes and compound dry–hot conditions across Iran, 1991–2024.”
- Rebuilt the manuscript around three explicit objectives; separated Results from Discussion; reduced repetition; added a 278-word quantitative abstract; calibrated causal, attribution, and long-return-period claims; and retained seven main figures and three evidence-focused tables.
- Used Alizadeh et al. (2020), supplied as `assets/Paper 3.pdf`, as a structural and rhetorical benchmark for empirical compound-event analysis. Its results-led organization and cautious event-component comparison informed the rewrite; no wording was copied.
- Corrected the threshold-window description to ±5 days (11 observations of calendar-day position before accounting for years) and added formal citations for Pettitt, standard normal homogeneity, and Buishand tests.
- Removed the exploratory composite fingerprint score from headline evidence because equal weights were not externally calibrated and several score components were dependent. Cluster assignments were demoted to supplementary, exploratory status because stability was uneven.
- Integrated the newly generated dependence-aware compound-trend results and internal temperature-consistency sensitivity into the Abstract, Methods, Results, Discussion, and Conclusions.
- Added `reports/Supplementary_Material_Q1.md`, which indexes 25 supporting tables and eight supporting figures and states the inferential limits of those outputs.
- Added `reports/Scientific_Audit_Q1.md`, documenting strengths, weaknesses, new analyses, claim calibration, reproducibility status, and a reviewer-style “major revision before submission” verdict.
- Added incomplete-but-visible submission declarations for data availability, code availability, author contributions, funding, competing interests, and acknowledgements. Seven bracketed placeholders remain because the required factual information is absent from the repository.
- Verified all seven main figure links and all 33 supplementary links, checked principal numerical claims against their CSV outputs, parsed `config.yaml`, compiled all Python modules, and ran `git diff --check` successfully.
- The complete historical pipeline was not rerun because the current Python environment lacks `statsmodels` and package installation was unavailable. The four targeted sensitivity outputs were executed against the existing derived data; a clean full rerun remains a pre-submission requirement.

## 2026-08-15 — Discussion and Conclusions expanded

- Expanded the Discussion in `reports/Manuscript_Q1_Revised.md` from 936 to approximately 3,150 words in response to the author's concern that it was too short for a Q1 article.
- Reorganized the Discussion into six analytical subsections covering the estimand and day-night asymmetry, regional heterogeneity, compound dry–hot interpretation, robustness and temporal uncertainty, methodological transferability, and research priorities.
- Added explicit distinctions between quantiles of annual event counts and quantiles of daily temperature, station-network extent and physical affected area, dry–hot conditions and drought, and observational concurrence and causal attribution.
- Expanded the Conclusions from 192 to approximately 540 words, separating the primary quantitative findings, spatial and climate-regime evidence, compound-event findings, robustness basis, and limits of inference.
- Applied a second-pass AI-writing audit using the `avoid-ai-writing` skill. Removed conversational and inflated transitions while retaining technical qualifiers required for scientific accuracy.

## 2026-08-15 — Reference-informed sentence and paragraph style revision

- Revised `reports/Manuscript_Q1_Revised.md` throughout to more closely match the rhetorical register of Alizadeh et al. (2020) while preserving scientific meaning, numerical values, citations, and claim strength.
- Merged short procedural and interpretive sentences where they represented a single logical unit, then restored selected short sentences to preserve the reference paper's alternation between compact statements and longer contrastive or causal constructions.
- Increased the controlled use of active scientific first person through constructions such as “we calculated,” “we estimated,” “our analysis shows,” and “we interpret,” without converting the manuscript into an author-centered narrative.
- Increased genuine causal, contrastive, and qualifying links using `because`, `while`, `whereas`, `although`, and `in contrast`; transition words were tied to explicit logical relations rather than added as decoration.
- Consolidated the Discussion from six subsections to four denser sections: distributional and regional expression, compound dry-hot conditions, robustness and methodological implications, and scope and research priorities. The resulting Discussion is approximately 2,480 words, retaining the expanded scientific content while improving paragraph continuity.
- After revision, approximate mean sentence lengths were 32 words in the Introduction, 33 in Methods, 27 in Results, 30 in Discussion, and 31 in Conclusions, compared with about 29, 31, 24, 28, and no separate conclusion, respectively, in the supplied reference article.
- Increased lexical similarity only at the level of shared disciplinary register and reporting verbs. An automated consecutive-word comparison found no distinctive sentence-level borrowing; the longest shared sequence contained six generic words.
- Reapplied the `avoid-ai-writing` audit, confirmed that prohibited promotional and template phrases were absent, compiled the Python modules, and passed `git diff --check`.

## 2026-08-15 — Reference-style sentence and paragraph revision

- Reworked `reports/Manuscript_Q1_Revised.md` to more closely match the rhetorical register and sentence rhythm of Alizadeh et al. (2020) while preserving numerical results, citations, and inferential limits.
- Merged sequences of short declarative sentences and introduced more explicit causal, contrastive, and qualifying structures using `because`, `whereas`, `although`, `while`, and `in contrast` only where the logical relationship was supported.
- Added limited scientific first-person constructions such as `we calculated`, `our analysis shows`, `we interpret`, and `we excluded`; passive and impersonal constructions were retained where the actor was not rhetorically important.
- Increased lexical similarity at the level of shared disciplinary register and reporting verbs, not distinctive wording. An automated consecutive-word comparison found no non-generic overlap longer than six words with the reference paper.
- Reorganized the Discussion from six subsections to four connected sections and reduced it from approximately 3,150 to 2,480 words by consolidating repeated qualifications rather than removing scientific content.
- After revision, approximate mean sentence lengths were 28.1 words in Results and 30.3 words in Discussion, compared with 23.7 and 27.7 words in the corresponding reference-paper sections. The revised text now alternates short evidentiary statements with longer interpretive sentences instead of using uniformly short prose.
- Re-ran the `avoid-ai-writing` audit, checked code syntax, and verified the repository diff for whitespace errors.
