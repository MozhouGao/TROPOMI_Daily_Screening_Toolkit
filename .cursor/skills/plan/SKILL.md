---
name: plan
description: >-
  Plan stage — route work into .specs/milestones or hand off to decide.
  Use when the user says plan, scope, milestone, or describes work before implementing.
disable-model-invocation: true
---

# plan

Route → clarify → write overview + issues → commit under `.specs/`.

```
plan                 # Auto-route
plan milestone       # Force milestone
```

Does **not** implement code. Does **not** mark plans ready — that is `plan-ready`.

## Phase 1: Route

| User intent | Path |
|-------------|------|
| Milestone / feature / fix / scoped delivery | This skill → `.specs/milestones/{slug}/` |
| Trade-off / ADR / “should we” | Hand off to `decide` |
| Approaches still wide open | Write a thin draft, then prefer `plan-explore` before locking |

Ask only if routing is truly ambiguous.

## Phase 2: Clarify

Ask only what blocks a useful plan. Prefer defaults:

- Goal in one paragraph
- In / out of scope
- Testable acceptance criteria
- Open questions / risks (explicit, not silent)

Do not interview for implementation trivia that `build` can discover.

## Phase 3: Write

Create (or update):

```
.specs/milestones/{slug}/overview.md   # from .specs/milestones/_template/overview.md
.specs/milestones/{slug}/issues.md     # from .specs/milestones/_template/issues.md
```

**overview.md must include:** goal, scope in/out, approach (cite ADRs by pointer),
testable acceptance criteria, risks/open questions table, references.

**issues.md:** trackable work items with status (`Todo` / `In Progress` / `Done` / `Blocked`).

Slug: kebab-case, stable, short.

## Phase 4: Commit

```bash
git add .specs/milestones/{slug}
git commit -m "docs: plan {slug}"
git push
```

Work on a worktree/branch — not `main` (constitution).

## Output contract

- [ ] Milestone folder exists with overview + issues
- [ ] Acceptance criteria are testable (usable by `build` / `tdd`)
- [ ] Open questions listed explicitly
- [ ] Committed (and pushed when remote exists)

## Next

| Situation | Next |
|-----------|------|
| Options / risks still fuzzy | `plan-explore` |
| Plan looks complete | `plan-ready` |
| Blocked on a decision | `decide`, then return |

## Scope discipline

Do not start `build`. Do not expand into unrelated milestones. Large design
choices become ADRs via `decide`, not essay dumps in the overview.
