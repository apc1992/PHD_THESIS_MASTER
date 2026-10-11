# SOURCE OF TRUTH & EVIDENCE CONTROL

Use the following hierarchy when sources disagree. A lower-ranked source may still be useful for identifying questions, but it must not override higher-ranked evidence without explicit justification.

1. **Validated experimental results**
2. **Final, QC'd computational outputs**
3. **Reproducible analysis/code that generates the outputs**
4. **Documented workflow/protocol descriptions**
5. **Current manuscript/thesis text**
6. **Informal notes / memory / planning documents**

## Important exceptions

- A manuscript can describe an analysis that was never actually executed; verify against code/output.
- A script can exist without being the script used for the final result; verify provenance.
- A figure can look plausible while containing an implementation or labeling error; trace it back to its source.
- Pending experiments are not evidence until results exist and are validated.

## Required traceability

For an important result, record as much of the following as is available:

**claim → sample/cohort → input → code/workflow → parameters → output → figure/table → final text**
