---
name: decide
description: >-
  Capture an architecture/product decision as an ADR under .specs/decisions/.
  Use when the user says decide, trade-off, should we, or ADR.
disable-model-invocation: true
---

# decide

Clarify options → write ADR → update index → link milestones → commit.

```
decide              # From conversation context
decide {slug}       # Named decision
```

## Phase 1: Frame

Capture:

- Decision question (one sentence)
- Forces / constraints
- Deadline or reversibility (cheap vs expensive to change later)

If this is really a delivery plan, not a choice among options → hand off to `plan`.

## Phase 2: Options

List 2–4 real options. For each: pros, cons, and what it optimizes for.
Avoid strawmen. “Do nothing” is allowed when it is a real option.

## Phase 3: Write ADR

Use `.specs/decisions/_template.md`:

```
.specs/decisions/{category}/{YYYY-MM-DD}-{slug}.md
```

Categories are lightweight folders you choose (`architecture`, `product`,
`process`, …). Keep paths stable once linked.

Status: `proposed` | `accepted` | `superseded`.

Prefer **accepted** only when the user (or clear prior agreement) chose.
Otherwise leave `proposed` and list what would make it accepted.

## Phase 4: Index + links

1. Update `.specs/decisions/_index.md` with id/title/status/path.
2. If a milestone depends on this decision, link it from that milestone’s
   overview References (and note blockers in open questions until accepted).

## Phase 5: Commit

```bash
git add .specs/decisions
git commit -m "docs: decide {slug}"
git push
```

## Output contract

- [ ] ADR file exists from the template
- [ ] Options and consequences are explicit
- [ ] `_index.md` updated
- [ ] Related milestones linked when applicable

## Next

| Status | Next |
|--------|------|
| `accepted` and a milestone was blocked on it | `plan` / `plan-ready` |
| `proposed` | User confirms, or gather evidence; do not pretend it is settled |
| Supersedes an older ADR | Mark old one `superseded` and point at the new path |

## Scope discipline

One decision per ADR. Do not implement the chosen option here — that is `build`.
Do not bury milestone task lists inside the ADR — that is `plan`.
