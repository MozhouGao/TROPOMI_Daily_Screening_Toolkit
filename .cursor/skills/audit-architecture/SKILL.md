---
name: audit-architecture
description: >-
  Advisory architecture and boundary audit. Use when reviewing placement,
  layering, or coupling. Does not merge or replace merge-pr architect.
disable-model-invocation: true
---

# audit-architecture

Advisory review of placement, layering, and coupling for a nominated scope
(PR diff, paths, or milestone). Not a substitute for `merge-pr`’s `architect`
on a shipping PR — use that crew when merging.

## Phases

1. Fix scope (diff / paths / milestone approach).
2. Map modules/packages and dependency directions relevant to the scope.
3. Check against the lens below.
4. Emit findings (do not auto-fix unless asked).

## Lens

- Change (or code) lives in the right layer/package
- Circular dependencies or layering violations
- Public surface larger than needed
- Parallel patterns that should share one approach
- Specs/ADRs missing when design has already drifted
- Premature platforms / abstractions (simplicity rule)

## Out of scope

| Concern | Where |
|---------|--------|
| Line-level lint/tests/debug | `audit-tests` / reviewer’s lane |
| Secrets / authz | security lane / `merge-pr` `security` |
| API schema compatibility detail | `audit-contracts` |

## Output

```text
id:       audit-architecture:{n}
severity: CRITICAL | HIGH | MEDIUM | LOW
file:     path
summary:  …
rule:     optional (layering, coupling, surface, follow-up, …)
fix:     direction only unless user asks to implement
```

End with a short “structure health” note: sound / frayed / needs replanning.
