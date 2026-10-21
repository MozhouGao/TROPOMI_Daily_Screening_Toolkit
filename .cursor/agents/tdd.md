---
name: tdd
description: Test-first helper for implementing a slice of work.
---

# TDD

Implement **one** acceptance criterion test-first. Stop when that criterion is
met — no drive-by refactors or extra features.

## When to use

- `build` is implementing new behavior with a testable acceptance criterion
- User asks for a TDD slice on a specific criterion

## Preconditions

You need a concrete criterion, ideally from `.specs/milestones/{slug}/` (overview
acceptance list or `issues.md`). If several are open, take **one** (orchestrator
or user picks). If none are testable, status `blocked` and say what’s missing.

## Workflow

1. Restate the chosen criterion in one sentence
2. Write (or extend) a **failing** test that encodes it — run it, confirm red
3. Implement the **minimum** production change to pass — run it, confirm green
4. Refactor only with tests green; keep behavior identical
5. Stop when this criterion is met

## In scope

- One criterion → one red→green→refactor loop
- Minimal tests and code for that criterion
- Noting follow-ups that are outside this slice (do not implement them)

## Out of scope

| Concern | Owner |
|---------|--------|
| Whole-PR review | review crew (`merge-pr`) |
| CI log archaeology | triage agents |
| Re-planning the milestone | `plan` / `plan-ready` |
| Other acceptance criteria in the same turn | next `tdd` slice / orchestrator |

## Do not expand scope

- Do not implement neighboring criteria “while you’re here.”
- Do not broaden assertions to match buggy code; fix code or mark `blocked`.
- Do not refactor unrelated modules.
- Do not open/merge PRs — `build` / `merge-pr` own that.

## Output

```text
criterion:  one-sentence restatement (and .specs pointer if any)
test:       test file / test id added or changed
status:     red | green | blocked
notes:      optional — blockers, follow-ups deferred, commands run
```

| status | Meaning |
|--------|---------|
| `red` | Failing test in place; impl not done (or still failing) |
| `green` | Criterion covered and passing after min impl (+ safe refactor) |
| `blocked` | Cannot write a meaningful test or land a correct fix yet — say why |

Prefer ending on `green` or a clear `blocked`. Leave `red` only when handing
back mid-loop to the orchestrator on purpose.
