---
name: audit
description: >-
  Advisory audit orchestrator. Use when the user asks for an audit without
  shipping code. Fans out to specialized audit skills as needed.
disable-model-invocation: true
---

# audit

Advisory only — does **not** merge, does **not** replace `merge-pr` review crew,
and does not change CI gate ownership.

```
audit                    # Route from user focus
audit architecture       # → audit-architecture
audit contracts          # → audit-contracts
audit tests              # → audit-tests
```

## Phase 1: Scope

Clarify:

- Target: whole repo, a path, a PR diff, or a milestone
- Focus: architecture / contracts / tests / combination
- Depth: quick pass vs thorough

Default target: current PR diff if one exists, else user-nominated paths.

## Phase 2: Route

| Focus | Skill |
|-------|-------|
| Architecture / boundaries / coupling | `audit-architecture` |
| Schemas / APIs / public contracts | `audit-contracts` |
| Test adequacy / regression locks | `audit-tests` |
| Unclear / “full audit” | Run all three (parallel when possible), then synthesize |

## Phase 3: Synthesize

Merge findings. Prefer the shared shape:

```text
id:       audit:{n}   # or audit-architecture:1 if single-skill run
severity: CRITICAL | HIGH | MEDIUM | LOW
file:     path
summary:  one or two sentences
rule:     optional tag
fix:     recommended direction (not a full implementation unless asked)
```

Order CRITICAL → LOW. Deduplicate overlaps (e.g. architecture vs contracts on the
same boundary).

## Rules

1. Do **not** auto-fix unless the user asks.
2. Do not expand into feature work or drive-by refactors.
3. Large remediations → follow-up issues / `plan`, not silent scope creep.
4. If the user wanted a merge gate, tell them to use `merge-pr` instead.

## Output

- Severity-ordered findings
- Short executive summary (what’s fine vs what’s not)
- Suggested next: fix now / file follow-ups / run `plan` / run `merge-pr`
