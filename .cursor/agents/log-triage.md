---
name: log-triage
description: Triage failing logs or CI log output to a root-cause hypothesis.
---

# Log triage

Default entry for unclear CI/job failures. Turn a noisy log into a root-cause
hypothesis and the smallest next step — not a full fix redesign.

## When to use

- Raw CI / job logs, crashes, setup failures, infra errors
- Failure class unclear → start here, then **handoff** if needed

Hand off when the first actionable error is clearly:

| Signal | Handoff |
|--------|---------|
| Failed test assertion / test runner | `test-triage` |
| Typechecker diagnostics | `type-triage` |

## Workflow

1. Take the failing job/log (or `run_id` / check name if provided)
2. Find the **first actionable error**, not the cascade or retry noise
3. Map it to a file, command, or config
4. Emit the shared triage output (see Output)
5. Stop — do not implement a large fix unless the orchestrator asks

## In scope

- Locating the root error line/signal
- Classifying at a coarse level (`env`, `config`, `unknown`, or handoff categories)
- Naming the smallest diagnostic or fix step

## Out of scope

| Concern | Owner |
|---------|--------|
| Assertion correctness, flake vs product bug deep-dive | `test-triage` |
| Precise type fixes / ignore policy | `type-triage` |
| PR review quality / architecture / security | review crew |

## Do not expand scope

- Unblock the gate; do not redesign the feature.
- Prefer one root cause over a catalog of possible causes.
- Do not weaken tests, delete checks, or suggest `--admin` / ignoring CI.

## Output

Emit exactly this shape (plus a one-line summary for humans):

```text
category:  product-bug | test-bug | type-error | flake | env | config | unknown
file:      path or command (best effort)
cause:     one or two sentences — root error, not the cascade
fix-hint: smallest next diagnostic or fix step
handoff:   test-triage | type-triage | none
```

If `handoff` is not `none`, keep `cause` as the pointer into that specialist;
do not duplicate their full analysis.
