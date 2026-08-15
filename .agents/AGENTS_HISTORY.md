\# AGENTS.md



\## Persistent Project History Protocol



This repository uses a persistent cumulative project memory file named



`History.md`



The purpose of `History.md` is to preserve the complete working history of the project across Codex sessions so that future sessions can reconstruct what happened, why decisions were made, what was changed, what failed, and what remains to be done.



Maintaining `History.md` is a mandatory part of every task.



\---



\# 1. Mandatory Startup Procedure



At the beginning of every Codex session, before making meaningful changes to the project



1\. Check whether `History.md` exists in the repository root.

2\. If it exists

&#x20;  - Read it before starting work.

&#x20;  - Recover relevant project context, previous decisions, unresolved issues, assumptions, completed tasks, and pending tasks.

&#x20;  - Treat it as persistent project memory.

3\. If it does not exist

&#x20;  - Create it using the structure defined below.

4\. Do not ask the user to repeat information already available in `History.md` unless that information is contradictory or genuinely ambiguous.



\---



\# 2. Core Rule History Must Be Cumulative



`History.md` is an append-only historical record.



Never



\- delete previous history;

\- overwrite the complete file;

\- truncate old entries;

\- silently rewrite previous decisions;

\- remove previous errors;

\- remove failed approaches;

\- erase abandoned experiments.



Previous entries represent the historical state at that moment and must remain preserved.



If an earlier decision later changes, append a new entry explaining



\- what changed;

\- why it changed;

\- which previous decision it supersedes.



Do not silently edit history to make the project appear cleaner than it actually was.



\---



\# 3. When History.md Must Be Updated



Update `History.md` whenever any meaningful interaction or project event occurs.



At minimum, update it after



\- every meaningful user request;

\- every meaningful clarification from the user;

\- every significant answer given to the user;

\- every implementation step;

\- every important investigation;

\- every architectural decision;

\- every change of approach;

\- every important assumption;

\- every file creation;

\- every meaningful file modification;

\- every file deletion or rename;

\- every important command;

\- every important tool call;

\- every test;

\- every validation;

\- every error;

\- every debugging attempt;

\- every discovered issue;

\- every resolved issue;

\- every important intermediate result;

\- every milestone;

\- every rollback;

\- every failed approach;

\- every important external dependency discovered;

\- every significant data-processing step;

\- every user correction;

\- every change requested by the user.



Do not wait until the entire task is complete.



Record important milestones immediately so that work can be recovered even if the session terminates unexpectedly.



\---



\# 4. Mandatory Pre-Response Synchronization



Before sending the final response for any meaningful task



1\. Review what happened during the current step.

2\. Update `History.md`.

3\. Make sure important decisions and modifications have been recorded.

4\. Record the current project state.

5\. Record unresolved issues and next steps.

6\. Only then provide the final response to the user.



The history synchronization is part of task completion.



A task is not considered complete until `History.md` has been updated.



\---



\# 5. Conversation Logging



For each meaningful user-agent exchange, preserve the interaction.



Record



\### User Request



Preserve the user's actual request as faithfully as practical.



For short or medium requests, record it verbatim.



For extremely long prompts, preserve



\- the complete intent;

\- explicit requirements;

\- constraints;

\- requested outputs;

\- important terminology;

\- corrections;

\- acceptance criteria.



Do not omit requirements merely to shorten the log.



\### Agent Response



Record the important content of the answer provided to the user.



When practical, preserve the final response verbatim.



For very large responses, record a detailed faithful summary containing all important



\- conclusions;

\- instructions;

\- decisions;

\- outputs;

\- caveats;

\- file references;

\- next actions.



Do not recursively copy the `History.md` logging instructions into `History.md`.



\---



\# 6. Reasoning Documentation



Private chain-of-thought, hidden reasoning, or internal scratchpad content must NOT be recorded.



Instead record a concise and useful Reasoning Summary.



The Reasoning Summary should explain



\- what problem was being solved;

\- what evidence was considered;

\- what assumptions were made;

\- what alternatives were considered when relevant;

\- why the selected approach was chosen;

\- why alternatives were rejected when relevant;

\- what uncertainty remains.



The objective is to make decisions auditable and reproducible without exposing private internal reasoning.



Example



```markdown

\### Reasoning Summary

The existing parser was retained because its output schema is already

used by three downstream modules. Replacing it would require unnecessary

changes to those modules. Input validation was therefore added around

the existing parser instead of replacing the parser itself.

```



Do not write meaningless statements such as



&#x20;I thought about the problem and decided this was best.



Document useful engineering reasoning.



\---



\# 7. History Entry Format



Every meaningful entry should use the following structure.



```markdown

\---



\## \[YYYY-MM-DD HHMMSS] — Short Event Title



Session session identifier if available



Status Planned  In Progress  Completed  Failed  Blocked  Revised



\### User Request



user request or faithful representation



\### Agent Response  Result



response, result, or detailed summary



\### Reasoning Summary



concise explanation of reasoning and decision basis



\### Decisions



\- decision

\- decision



\### Assumptions



\- assumption

\- assumption



\### Actions Performed



\- action

\- action



\### Files Changed



\- `pathtofile` — created  modified  renamed  deleted

\- `pathtoanotherfile` — modified



\### Commands  Tools



\- `important command`

&#x20; - Exit status `code if available`

&#x20; - Result short useful result



\### Tests  Validation



\- test

\- Result PASS  FAIL  PARTIAL  NOT RUN



\### Errors  Problems



\- error or problem

\- Cause known or suspected cause

\- Resolution resolution if available



\### Current State



state of the relevant part of the project after this step



\### Next Steps



1\. next action

2\. next action



\---

```



Sections with no meaningful information may be omitted.



Do not fill sections with useless text such as `NA` unless the absence itself matters.



\---



\# 8. Session Tracking



At the beginning of a new working session append a session marker



```markdown

\---



\# Session Started — YYYY-MM-DD HHMMSS



\## Recovered Context



\- Current project objective

\- Previous completed work

\- Important existing decisions

\- Current unresolved issues

\- Immediate next step



\---

```



If a reliable session identifier is available, include it.



Otherwise use a timestamp-based identifier such as



`session-2026-08-13-095300`



\---



\# 9. Current Project State



At important milestones, append a checkpoint describing the current state.



Use



```markdown

\## Project Checkpoint — YYYY-MM-DD HHMMSS



\### Completed



\- ...



\### In Progress



\- ...



\### Pending



\- ...



\### Known Problems



\- ...



\### Important Decisions



\- ...



\### Important Files



\- `...`

\- `...`



\### Recommended Next Action



...

```



This makes recovery easier without replacing the detailed chronological history.



\---



\# 10. User Corrections Have High Priority



If the user corrects



\- a requirement;

\- an assumption;

\- a parameter;

\- a methodology;

\- a filename;

\- a path;

\- a design decision;

\- a previous interpretation;



record the correction explicitly.



Example



```markdown

\### User Correction



Previous interpretation

Use daily aggregation.



Corrected requirement

Use hourly data without daily aggregation.



Impact

The preprocessing and validation pipeline must retain hourly resolution.

```



Future work must follow the corrected requirement.



\---



\# 11. Decision Changes



Never hide changed decisions.



When changing a previous decision, record



```markdown

\### Decision Revision



Previous decision

...



New decision

...



Reason for revision

...



Affected filescomponents

...

```



This is especially important for scientific, engineering, analytical, and software-development projects.



\---



\# 12. File Change Tracking



For every meaningful file modification record



\- file path;

\- action;

\- purpose.



Example



```markdown

\### Files Changed



\- `srcpreprocess.py` — modified

&#x20; - Added missing-value validation.

&#x20; - Added CRS consistency check.



\- `teststest\_preprocess.py` — created

&#x20; - Added tests for missing values and invalid CRS.

```



Do not paste entire files into `History.md` unless explicitly necessary.



History should describe the change, while Git or the working tree contains the actual implementation.



\---



\# 13. Command Tracking



Record commands when they materially affect the project or help reproduce a result.



Examples



```markdown

\### Commands  Tools



\- `pytest teststest\_preprocess.py -q`

&#x20; - Exit status 0

&#x20; - Result 14 tests passed.



\- `python scriptsdownload\_data.py`

&#x20; - Exit status 1

&#x20; - Result Download failed because API credentials were missing.

```



Do not record every trivial shell command.



Commands such as



```bash

ls

pwd

cat

```



only need to be recorded when their result materially influences a decision.



\---



\# 14. Error and Debugging History



Failures are important project knowledge.



Never remove them from history.



For meaningful errors record



```markdown

\### Error



Observed

error



Context

what operation caused it



Likely Cause

cause



Attempted Fixes

1\. ...

2\. ...



Resolution

...



Status

Resolved  Unresolved  Workaround

```



Failed approaches must remain documented so they are not repeated blindly in future sessions.



\---



\# 15. Scientific and Analytical Work



For scientific, research, statistical, modeling, GIS, climate, engineering, or data-analysis tasks, additionally record whenever relevant



\- dataset;

\- dataset source;

\- spatial resolution;

\- temporal resolution;

\- study period;

\- variables;

\- units;

\- coordinate system;

\- preprocessing;

\- quality-control operations;

\- missing-data handling;

\- equations;

\- parameters;

\- parameter sources;

\- assumptions;

\- thresholds;

\- models;

\- software;

\- package versions when important;

\- statistical tests;

\- validation methods;

\- evaluation metrics;

\- output figures;

\- output tables;

\- output maps;

\- uncertainty;

\- limitations.



Any parameter that materially affects results should have its origin recorded.



Example



```markdown

\### Parameter Provenance



\- `threshold = 0.38`

&#x20; - Purpose glare acceptability threshold

&#x20; - Source project methodology  referenced standard

&#x20; - User-provided No

&#x20; - Hard-coded No

```



Do not allow unexplained magic numbers to silently enter the project.



\---



\# 16. Data Provenance



When data is introduced into the project, record



```markdown

\### Data Provenance



Dataset

...



Source

...



URL  DOI  repository

...



Acquisition method

Manual  API  script  existing local file



Original files

...



Processing steps

...



Generated files

...



Notes

...

```



If data was provided directly by the user, state that explicitly.



\---



\# 17. Testing and Validation



Important code or analysis changes must include the validation status.



Use one of



\- `PASS`

\- `FAIL`

\- `PARTIAL`

\- `NOT RUN`



Example



```markdown

\### Validation



Unit tests PASS



Integration tests PASS



Scientific result validation PARTIAL



Reason

Implementation was verified against the reference equations, but no

independent observational dataset is currently available.

```



Never imply validation occurred if it did not.



\---



\# 18. Uncertainty and Unverified Claims



If something has not been verified, explicitly record it.



Use phrases such as



\- `Unverified`

\- `Assumed`

\- `Requires confirmation`

\- `Inferred from available evidence`



Never turn assumptions into historical facts.



\---



\# 19. Security and Sensitive Information



Never store secrets in `History.md`.



Redact



\- passwords;

\- API keys;

\- authentication tokens;

\- private keys;

\- cookies;

\- credentials;

\- access tokens;

\- connection secrets;

\- confidential environment variables.



Use



`\[REDACTED]`



Example



```markdown

API\_TOKEN=\[REDACTED]

```



Do not copy secret-containing command output into history.



\---



\# 20. Logging Noise Control



`History.md` should be comprehensive but useful.



Do not log



\- trivial navigation;

\- repetitive directory listings;

\- insignificant formatting changes;

\- repeated unchanged observations;

\- every token of terminal output;

\- enormous dependency installation logs;

\- generated binary contents;

\- hidden chain-of-thought.



Summarize large command outputs while preserving



\- command;

\- exit status;

\- important result;

\- important error messages;

\- relevant paths;

\- key metrics.



\---



\# 21. Large Output Handling



If command output is large, record a concise summary.



Example



Instead of storing 5,000 lines of test output



```markdown

\### Tests



Command

`pytest -q`



Result

PASS



Summary

187 tests passed in the full test suite.



Warnings

3 deprecation warnings from dependency X.

```



If a specific error line matters, preserve that line.



\---



\# 22. Git Awareness



When Git is available, History.md complements Git but does not replace it.



When relevant record



\- branch;

\- commit hash;

\- changed files;

\- meaningful diff summary.



Do not assume a commit exists unless it actually exists.



Never create a commit solely because this protocol exists unless the user requested commits or repository workflow requires them.



\---



\# 23. Before Risky Operations



Before performing an operation that could significantly modify or destroy project state, create a checkpoint in `History.md`.



Examples



\- major refactoring;

\- mass file modification;

\- database migration;

\- data deletion;

\- dependency replacement;

\- architectural rewrite;

\- generated-output replacement.



Record



```markdown

\## Pre-Change Checkpoint



Planned operation

...



Reason

...



Affected components

...



Expected result

...



Rollback strategy

...

```



\---



\# 24. Automatic Recovery Behavior



If a new session starts with little conversational context



1\. Read `History.md`.

2\. Inspect relevant project files.

3\. Recover the latest confirmed project state.

4\. Continue from that state.



Do not restart the project from scratch unless explicitly requested.



When there is a conflict between



\- current user instruction;

\- `History.md`;

\- current repository state;



use this priority



1\. Current explicit user instruction

2\. Current repository state when objectively verifiable

3\. Latest confirmed entry in `History.md`

4\. Older historical entries



Record the conflict and its resolution.



\---



\# 25. Checkpoint Command



If the user sends



`checkpoint`



immediately append a comprehensive project checkpoint to `History.md`.



It must include



\- what has been done;

\- current state;

\- files changed;

\- major decisions;

\- unresolved errors;

\- assumptions;

\- validation status;

\- pending tasks;

\- exact recommended next step.



Do not perform unrelated work for a `checkpoint` request unless the user asks for it.



\---



\# 26. History Requests



If the user asks questions such as



\- What have we done

\- Where did we stop

\- Why did we do this

\- What changed

\- What errors did we have

\- What files did you modify

\- What remains

\- What did I ask you before



consult `History.md` before answering.



Prefer documented history over guessing from memory.



\---



\# 27. History.md Initial Template



When `History.md` does not exist, initialize it with



```markdown

\# Project History



&#x20;Persistent cumulative development and research history.



&#x20;This file is maintained automatically according to `AGENTS.md`.

&#x20;Historical entries must not be deleted or silently rewritten.



\## Project



Repository repository name if known



History Started timestamp



Purpose project purpose if known



\---



\## History Conventions



\- Entries are chronological.

\- Previous entries are preserved.

\- Corrections are recorded as new entries.

\- Sensitive credentials are redacted.

\- Reasoning is recorded as concise decision summaries, not private chain-of-thought.



\---



\# Session Started — timestamp



\## Initial Context



available project context



\---

```



\---



\# 28. End-of-Task Entry



Before finishing a significant task, append



```markdown

\## Task Completion — YYYY-MM-DD HHMMSS



\### Requested



...



\### Delivered



...



\### Files Changed



\- ...



\### Validation



...



\### Remaining Issues



\- ...



\### Current State



...



\### Recommended Next Step



...

```



\---



\# 29. Never Falsify History



Never claim



\- a command was executed when it was not;

\- a test passed when it was not run;

\- a file was modified when it was not;

\- a dataset was validated when it was not;

\- a source was consulted when it was not;

\- an error was resolved when only a workaround exists.



`History.md` is intended to be an auditable record.



Accuracy is more important than making the project appear successful.



\---



\# 30. Mandatory Operating Loop



For every meaningful task follow this loop



```text

READ History.md

&#x20;       ↓

UNDERSTAND current project state

&#x20;       ↓

RECORD user request

&#x20;       ↓

INVESTIGATE  IMPLEMENT

&#x20;       ↓

RECORD important intermediate events

&#x20;       ↓

TEST  VALIDATE

&#x20;       ↓

RECORD results and decisions

&#x20;       ↓

UPDATE current state + next steps

&#x20;       ↓

SYNC History.md

&#x20;       ↓

RESPOND TO USER

```



For multi-step tasks, repeat the recording process at meaningful milestones rather than waiting until the end.



\---



\# 31. Highest-Priority History Rules



If there is any ambiguity about this protocol, follow these rules



1\. `History.md` must survive across sessions.

2\. Read it at the start of every session.

3\. Never erase historical entries.

4\. Record every meaningful user request.

5\. Record every meaningful result.

6\. Record important decisions and their rationale.

7\. Record file changes.

8\. Record important commands and tests.

9\. Record failures as well as successes.

10\. Redact secrets.

11\. Use concise reasoning summaries instead of private chain-of-thought.

12\. Update `History.md` before the final response.

13\. Leave enough information that another Codex session can continue the project without requiring the user to explain everything again.



These requirements apply throughout the repository unless a more specific nested `AGENTS.md` explicitly overrides a particular rule.

