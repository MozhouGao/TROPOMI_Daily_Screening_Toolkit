---
name: build
description: >-
  Implement on a draft PR until local gates are green. Use when the user says
  build, implement, or continue a milestone slice after plan-ready or scaffold.
  Dispatches tdd per testable acceptance criterion; runs verify; hands off to
  merge-pr. Leaves PRs as draft — does not run the review crew.
disable-model-invocation: true
---

# build

Implement behavior on an issue branch / worktree until local verification is
green. Open or update a **draft** PR; do **not** mark it ready or merge — that
is `merge-pr`.

```
build                    # Current branch / milestone slice
build {milestone-slug}   # Named milestone under .specs/milestones/
build --criterion "…"    # One acceptance criterion only (single tdd slice)
```

**Ownership:** `build` implements and verifies locally. `merge-pr` is the sole
review-crew dispatch point. Do not run `reviewer` / `architect` / `security` /
`adversary` here.

---

## Phase 1: Preconditions

| Check | Fail action |
|-------|-------------|
| Not on `main` / `master` | Stop; create worktree + issue branch per `AGENTS.md` |
| Plan exists (`.specs/milestones/{slug}/`) or user gave explicit scope | If missing → `plan` or `plan-ready`, not guesswork |
| Prefer **GO** from `plan-ready` (or explicit user waiver) | Note waiver in PR body if skipping |
| Working tree understood (`git status`) | Commit or stash unrelated WIP before starting a slice |

Read the milestone overview, issues, and acceptance criteria. List open,
testable criteria for this session.

---

## Phase 2: Branch and draft PR

```bash
git fetch origin
git status --porcelain
gh pr view --json number,title,state,isDraft,url 2>/dev/null || true
```

| Situation | Action |
|-----------|--------|
| No open PR for branch | Push branch; `gh pr create --draft` using `.github/pull_request_template.md` |
| Open draft PR exists | Continue on it; keep `isDraft: true` |
| Open **ready** PR | Leave ready state alone unless user asks; prefer a new branch for new scope |
| No remote branch | `git push -u origin HEAD` before creating PR |

PR body should link the milestone path and list criteria in scope / done / deferred.

---

## Phase 3: TDD slices (sole `tdd` dispatch point)

For each testable acceptance criterion (one per slice unless user narrows scope):

1. **Load** `.cursor/agents/tdd.md` and run with that file as instructions plus:
   - The chosen criterion (from `.specs/` or user)
   - Relevant stub paths from `scaffold` if any
   - Prior slice outputs in this session
2. Collect the `tdd` output contract:
   `criterion`, `test`, `status` (`red` \| `green` \| `blocked`), `notes`
3. On `green` → commit meaningful progress (Conventional Commits, ≤ 120 chars)
4. On `blocked` → stop the build loop; report blocker and suggested next skill
   (`plan`, `decide`, `scaffold`)
5. On `red` mid-loop → either continue the slice or hand back to orchestrator;
   do not start the next criterion until the current one is `green` or `blocked`

Dispatch slices **sequentially** — one criterion per `tdd` invocation.

### Scope discipline

- One criterion per `tdd` slice; no drive-by refactors or extra features
- Match existing repo conventions (see workspace rules)
- New behavior ships with tests (workspace testing rule)
- Update `.specs/` on the same branch when discoveries affect the plan (`docs:` commits)

---

## Phase 4: Verify (local gates)

Run the `verify` skill after all planned slices are `green` (or when the user
asks to check gates mid-build).

```bash
# verify skill runs the configured gate table (lint, typecheck, test, …)
```

| `verify` overall | Action |
|------------------|--------|
| `GREEN` | Continue to Phase 5 |
| `RED` | Fix via triage agents (see `verify` → **Triage routing**) or code; re-run `verify` |
| `BLOCKED` | Stop; report environment blocker; do not claim done |

`build` does not replace CI — `merge-pr` still requires green PR checks.

---

## Phase 5: Push and handoff

```bash
git status --porcelain
git push
gh pr view --json number,isDraft,url
```

**Done when:**

- [ ] In-scope acceptance criteria are `green` (or explicitly deferred in PR body)
- [ ] Changes committed and pushed
- [ ] Draft PR open and up to date
- [ ] `verify` overall is `GREEN`
- [ ] `.specs/` updated if the implementation changed the plan
- [ ] PR remains **draft**

### Output template

```markdown
**Build**: {branch} → draft PR #{n} ({url})
**Milestone**: `.specs/milestones/{slug}/`
**Slices**: {n} tdd — {green count} green, {blocked count} blocked
**Verify**: GREEN | RED | BLOCKED
**Deferred**: {criteria or none}

**Next**: `merge-pr` when ready to ship (review crew + CI + merge)
```

---

## What build does not do

- Run the review crew (`merge-pr` only)
- Mark PR ready or merge
- Expand milestone scope without `plan` / user approval
- Skip tests for new behavior
- Implement on `main`

---

## Mid-milestone resume

When continuing a milestone across sessions:

1. Read the draft PR body and `.specs/milestones/{slug}/`
2. Identify remaining open criteria
3. Resume at Phase 3 for the next criterion
