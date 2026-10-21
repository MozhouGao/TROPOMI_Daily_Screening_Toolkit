---
name: test-triage
description: Triage failing tests to root cause vs flake vs bad assertion.
---

# Test triage

Triage a failing test into a clear category and the smallest correct fix path.
Prefer truth over making the suite green.

## When to use

- Failed assertions / test-runner output (pytest, jest, go test, etc.)
- Handed off from `log-triage` when the first error is clearly a test failure

If the log is mostly infra/setup and never reaches an assertion → hand back
toward `log-triage` (`handoff` is not used for that; set `category: env` or
`config` and `fix-hint` accordingly).

## Workflow

1. Identify the failing test id/name and assertion
2. Relate it to the PR diff (or state that it looks unrelated/env)
3. Classify (see Categories)
4. Recommend a fix path that does **not** casually weaken assertions
5. Emit the shared triage output (see Output)

## Categories

| category | Meaning |
|----------|---------|
| `product-bug` | Code under test is wrong |
| `test-bug` | Test is wrong, outdated, or over-specified vs intended behavior |
| `flake` | Non-deterministic; see Flake rule |
| `env` | Tooling, fixtures, services, credentials, platform |
| `config` | Misconfigured test harness / paths / env flags |
| `unknown` | Cannot classify yet — say what signal is missing |

## Flake rule

Call `flake` only with evidence. Align with `merge-pr`:

- One fail then one pass on the **same SHA** → may note possible flake; still investigate
- Prefer **two failures on the same SHA** before treating flake as the working theory
- Never “fix” flake by deleting the test or broadening assertions without cause

## In scope

- Root cause of the failing test
- Whether to fix product code vs test vs environment
- Minimal fix-hint

## Out of scope

| Concern | Owner |
|---------|--------|
| Raw job setup / infra before tests run | `log-triage` |
| Typechecker-only failures | `type-triage` |
| Broad PR review | review crew |

## Do not expand scope

- Do not rewrite the feature to dodge a single assertion.
- Do not recommend weakening or removing assertions unless `category: test-bug` and the assertion is demonstrably wrong.
- Do not expand into unrelated test cleanup.

## Output

```text
category:  product-bug | test-bug | type-error | flake | env | config | unknown
file:      test file and/or code-under-test path
cause:     failing assertion + root hypothesis (one or two sentences)
fix-hint: smallest correct fix path
handoff:   type-triage | none
```

Use `handoff: type-triage` only when the “test” failure is actually a typecheck
step mislabeled in CI. Otherwise `none`.
