---
name: audit-tests
description: >-
  Advisory test-adequacy audit. Use when reviewing coverage of new behavior,
  flaky tests, or missing regression locks. Does not merge or replace merge-pr.
disable-model-invocation: true
---

# audit-tests

Advisory review: do tests lock acceptance criteria and important regressions?

## Phases

1. Fix scope (PR diff, paths, or milestone acceptance list).
2. Map criteria / changed behavior → existing tests.
3. Check against the lens below.
4. Emit findings (do not auto-fix unless asked).

## Lens

- Each testable acceptance criterion has (or clearly lacks) a lock
- Regressions adjacent to the change are covered when risk is high
- Assertions are meaningful (not tautologies / over-broad mocks)
- Flaky patterns: time, network, order dependence, shared mutable state
- Fixtures are synthetic / safe (no secrets or prod data)
- Missing cases listed before recommending new tooling or coverage % targets

## Out of scope

| Concern | Where |
|---------|--------|
| Package layering | `audit-architecture` |
| Public contract shape / versioning | `audit-contracts` |
| Implementing the missing tests | `tdd` / `build` (only if user asks to fix) |

## Output

```text
id:       audit-tests:{n}
severity: CRITICAL | HIGH | MEDIUM | LOW
file:     test or code path
summary:  …
rule:     optional (missing-lock, flake-risk, weak-assertion, follow-up, …)
fix:     direction only unless user asks to implement
```

Prefer “missing case X” over “increase coverage to N%” unless the repo already
enforces a coverage gate.
