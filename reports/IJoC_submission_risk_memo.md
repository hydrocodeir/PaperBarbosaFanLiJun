# IJoC Submission Risk Memo

## Overall judgement

The manuscript is now substantially safer than it was in earlier drafting stages and is meaningfully closer to a defensible `International Journal of Climatology` submission. Its main strengths are now clear: a focused central claim, a coherent quantile-regression framework, stronger discipline in separating main-text figures from supplementary diagnostics, and a more explicit robustness architecture. It no longer reads like an exploratory regional report with many outputs; it reads much more like a hypothesis-driven climatology manuscript with a structured analytical workflow.

The remaining risks are therefore narrower and more technical. They now concern the depth of selected methodological choices rather than the overall scientific seriousness of the paper. In reviewer terms, the manuscript is no longer most vulnerable to "Why is this paper built this way?" but to a smaller number of sharper questions such as "How should the homogeneity diagnostics be interpreted?", "Are some physical interpretations still stronger than the direct evidence?", and "Are any communication choices still open to over-reading?"

## Updated risk ratings

### 1. Bootstrap depth in the tails
**Risk level:** Low to moderate

This risk has now been reduced further. In addition to the existing methodological justification, the manuscript now includes a targeted higher-replication sensitivity check for the two central warm-extreme indices, increasing the bootstrap depth from `200` to `400` replicates for the core warm-day and warm-night summaries. The result is reassuring: regional bootstrap means changed only minimally, and station-level correspondence between the 200- and 400-replicate summaries remained extremely high.

**Residual concern**
The residual vulnerability is no longer "no evidence exists", but rather that the higher-depth check is targeted rather than universal. A very demanding reviewer could still ask why the same exercise was not extended to all four indices or to even deeper replication in the far tails.

**Best mitigation**
No additional work is essential. If more effort is available, the strongest further upgrade would be to cite the targeted high-replication result prominently in any response letter or reviewer response.

### 2. Data homogeneity and station quality
**Risk level:** Moderate

This risk has been materially reduced. The manuscript now documents screening for completeness, duplicates, and physical consistency; it explains the use of detrended Pettitt-, SNHT-, and Buishand-style diagnostics; and it includes an explicit exclusion sensitivity analysis. That is a major improvement. The manuscript can now show that the principal directional conclusions survive exclusion of flagged stations.

**Residual concern**
The residual vulnerability is interpretive rather than structural. A reviewer may still ask why annual mean temperature was used as the screening series, why flagged stations were retained in the main analysis, or whether a full homogenization workflow would materially alter the findings. The sensitivity analysis reduces this risk, but it also reveals a real nuance: warm-night tail asymmetry is not completely insensitive to record screening.

**Best mitigation**
No major new analysis is essential. The most useful further improvement would be a short supplementary note or one additional sentence in the manuscript clarifying why annual mean temperature was used as the common screening target and why the sensitivity result is sufficient for the manuscript's inferential scope.

### 3. Physical interpretation remains somewhat lighter than the statistical architecture
**Risk level:** Moderate

This has improved because the manuscript now more clearly distinguishes supported inference from plausible mechanism. The wording is more disciplined and less causal than before. Even so, the discussion still leans more on climatological plausibility than on direct process attribution, particularly for nocturnal heat retention, elevation-related contrasts, and inland-versus-maritime modulation.

**Why this is still a reviewer issue**
`IJoC` reviewers often tolerate regional observational papers with limited process attribution, but they are less tolerant when mechanism language sounds firmer than the evidence warrants.

**Best mitigation**
Mainly editorial. The remaining task is to ensure that no sentence anywhere in the paper drifts back into overconfident causal wording.

### 4. The slope-translated "approximate 1961-2024 change" column may still attract scrutiny
**Risk level:** Moderate

This item is better handled than before because the caption now clarifies that the quantity is a slope-based translation rather than an endpoint-to-endpoint observed difference. Even so, some reviewers dislike such columns because they can be over-read or mistaken for realized changes in a nonstationary record.

**Best mitigation**
If you want maximum conservatism, move that column to the Supplementary Material. If you keep it in the main text, the current wording is probably adequate but not risk-free.

### 5. Regionalization and clustering justification
**Risk level:** Low

This risk is now lower still. The manuscript not only explains why hierarchical clustering was chosen and how uncertainty-aware features were used, but now also includes an explicit alternative-method sensitivity analysis comparing the baseline partition with complete-linkage, Ward-linkage, cityblock-distance, and k-means alternatives. That additional evidence substantially reduces the risk that a reviewer would see the regionalization as depending on one arbitrary clustering convention.

**Residual concern**
A reviewer could still ask whether cooling-related regionalization is less stable than warming-related regionalization, but that is now something the manuscript can answer directly rather than avoid.

**Best mitigation**
Already addressed. No further methodological work is necessary here.

### 6. Split-period analysis
**Risk level:** Low

This risk has been reduced effectively. The manuscript now states explicitly that the `1961-1990` versus `1991-2024` comparison is a synthesis-oriented nonlinearity diagnostic rather than a formal breakpoint attribution exercise. That framing is methodologically honest and should diffuse the most obvious reviewer concern.

**Residual concern**
A reviewer may still ask why `1990` was chosen, but this is now unlikely to become a major obstacle unless the section is overstated elsewhere.

**Best mitigation**
No new analysis is necessary.

### 7. Figure density and manuscript presentation
**Risk level:** Low

This issue was previously more serious. It is now much better controlled because structurally rich diagnostics have been moved to the Supplementary Material and the main text has a clearer figure strategy. The remaining risk is mostly editorial: captions must stay lean, and the main text should not repeat every message already visible in the figures.

**Best mitigation**
A final pre-submission tightening pass on captions and local transitions is sufficient.

### 8. Citation hygiene and submission polish
**Risk level:** Low to moderate

The internal citation structure is much cleaner than before, and ambiguities such as `2005a/2005b` have been resolved. What remains here is ordinary submission polish: exact journal style, punctuation consistency, and ensuring that every interpretive claim has a citation where a reviewer expects one.

**Best mitigation**
Do one final journal-style cleanup immediately before submission.

## Most likely reviewer comments now

If the manuscript were submitted in its current form, the most plausible substantive reviewer comments are now likely to be:

1. "Please clarify why annual mean temperature was used for homogeneity screening and how the retained flagged stations affect warm-night tail interpretation."
2. "Please keep the physical discussion clearly inferential and avoid implying direct process attribution where the evidence is indirect."
3. "Please clarify again that the approximate full-period change is a slope-based translation rather than an observed endpoint difference."
4. "Please explain why the higher-replication bootstrap check was targeted to the warm indices and not extended to every analysis component."
5. "Please clarify how the more method-sensitive cool-day clustering should be interpreted."

These are all substantially more manageable than the earlier structural criticisms the manuscript would have attracted before the recent revisions.

## Probability judgement

If submitted now:

1. Probability of desk rejection: low
Reason: the paper now has a clear scientific question, a nontrivial analytical design, and a presentation style that is much closer to a real journal submission.

2. Probability of major revision: moderate
Reason: reviewers are still likely to push on homogeneity interpretation, the balance between statistical evidence and physical explanation, and a few communication choices in the main tables.

3. Probability of minor revision directly: low to moderate
Reason: possible if reviewers are sympathetic to the manuscript's observational-regional framing and do not insist on deeper tail uncertainty analysis.

4. Probability of eventual acceptance after revision: good to very good
Reason: the core scientific signal is coherent, the paper is now much more professionally framed, and the remaining weaknesses are concentrated rather than fundamental.

## Best next actions before submission

If more work is possible, the highest-yield actions are now:

1. Decide whether the slope-translated `Approx. 1961-2024 change` column should remain in `Table 1` or move to the Supplementary Material.
2. Perform one more line-edit pass focused only on causal phrasing so that all mechanism language remains explicitly inferential.
3. Add one short clarifying sentence, either in the main text or supplementary material, explaining why annual mean temperature was used as the common homogeneity-screening series.
4. Do a final journal-style cleanup of references, captions, metadata, and supplementary cross-references.
5. Prepare a concise reviewer-response note highlighting the new homogeneity-exclusion, higher-replication bootstrap, and alternative-clustering sensitivities.

## Bottom line

The manuscript is now past the stage where its main risks arise from structural weakness. The largest previously exposed vulnerability, namely the absence of an explicit homogeneity-exclusion sensitivity check, has already been reduced materially. Likewise, clustering justification and split-period framing are now much better defended, and the targeted higher-replication bootstrap analysis has reduced what had been the clearest remaining technical criticism. At this point, the paper's remaining exposure is concentrated less in missing analyses than in a small set of interpretive and presentation choices. That is a much healthier submission position than before.

