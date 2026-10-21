---
name: merge-pr
description: >-
  Ship stage — full review-and-merge pipeline for the current (or named) PR.
  Use when the user says merge-pr, merge this PR, ship it, review and merge,
  check CI and merge, or after build leaves a green draft. Sole review-crew
  dispatch point: Wave 1 reviewer + architect + security (when relevant);
  Wave 2 adversary; severity triage; CI green; draft→ready; squash-merge;
  cleanup. Never merges with failing checks; never uses --admin.
disable-model-invocation: true
---

# merge-pr

Merge a PR through the full review-and-merge pipeline. Accepts a PR number or
operates on the current branch's open PR.

```
merge-pr              # Current branch's PR
merge-pr 1234         # Specific PR number
merge-pr --dry-run    # Reviews + triage only; do not merge
```

**Ownership:** `build` leaves PRs as **draft** on purpose. This skill is the
canonical (and only) place that dispatches the review crew and flips
draft → ready → merge.

---

## Phase 1: Identify the PR

```bash
gh pr view --json number,title,state,headRefName,baseRefName,isDraft,mergeable,url
```

If given a number: `gh pr view {PR#} --json …`.

**Gate:** if there is no open PR for the branch, stop and say so. Do not invent one.

---

## Phase 2: Pre-flight checks

```bash
git status --porcelain
git fetch origin
git log --oneline origin/{base}..HEAD
gh pr view {PR#} --json mergeable,mergeStateStatus,headRefOid,baseRefName
```

Use the PR's `baseRefName` (usually `main`) as `{base}`.

**All must pass before Phase 3:**

| Check | Fail action |
|-------|-------------|
| Working tree clean | Stop; ask user to commit or stash |
| ≥ 1 commit ahead of base | Stop; PR is empty |
| `mergeable` is `MERGEABLE` | If `UNKNOWN`, re-poll once after a few seconds. If still not mergeable, stop |
| Not `CONFLICTING` / `DIRTY` | Stop; user must rebase/merge base and resolve. Do **not** dispatch the crew — reviewers would see the wrong diff |
| `BEHIND` | Allowed, but note it in a PR comment; prefer updating from base before merge |

---

## Phase 3: Review crew (sole dispatch point)

`build` does **not** run the review crew. Always dispatch here on the current
head SHA. Load persona prompts from `.cursor/agents/{name}.md` and run each
agent with that file as its instructions plus the PR diff context.

### Wave 1 — parallel

| Agent | Prompt | Mandate |
|-------|--------|---------|
| **reviewer** | `.cursor/agents/reviewer.md` | Code quality on the diff |
| **architect** | `.cursor/agents/architect.md` | Placement, boundaries, coupling |
| **security** | `.cursor/agents/security.md` | Secrets / trust boundaries — **only if** the diff touches auth, credentials, crypto, network exposure, secret-bearing config, or untrusted input at a boundary. Skip for pure docs/specs and note the skip |

Dispatch Wave 1 in parallel when the host supports it.

**Finding IDs:** when collecting Wave 1 output, assign each finding a stable id
(`reviewer:1`, `architect:2`, …) so Wave 2 can cross-reference them.

### Wave 2 — after Wave 1 completes

| Agent | Prompt | Mandate |
|-------|--------|---------|
| **adversary** | `.cursor/agents/adversary.md` | Cross-check Wave 1; regressions, symptom-only fixes, contradictions. Runs **last** |

Pass the id'd Wave 1 findings into the adversary prompt.

### Reviewer output contract

Each agent returns:

- **Verdict:** `APPROVE` | `APPROVE_WITH_CHANGES` | `BLOCK` | `SKIP`
- **Findings:** `{ id, severity, file, line?, summary, rule? }`
- Severities: `CRITICAL` > `HIGH` > `MEDIUM` > `LOW`

---

## Phase 3.5: Capture crew verdict (best-effort)

After Wave 2 and **before** triage fixes, snapshot crew output as emitted into
one PR comment. Capture must **never** block a merge — on failure, warn and continue.

Use a sentinel so the comment can be upserted:

```html
<!-- merge-pr:review-crew v1 -->
…json or compact markdown summary…
<!-- /merge-pr:review-crew -->
```

Include `headSha` and UTC `ranAt` from:

```bash
gh pr view {PR#} --json headRefOid --jq .headRefOid
```

Do not fabricate agent verdicts.

---

## Phase 4: Triage findings

| Severity | Action |
|----------|--------|
| **CRITICAL** | Fix before continuing |
| **HIGH** | Fix before continuing |
| **MEDIUM** | Fix if small; otherwise note follow-up in the PR body |
| **LOW** | Comment only unless trivial |

After fixes: commit, push. Re-run **only** agents whose findings were addressed
— not the full crew — unless the fix substantially changes the diff surface.

Any agent `BLOCK` verdict with unresolved CRITICAL/HIGH → do not proceed to CI merge.

---

## Phase 5: CI gate

No discretion.

```bash
gh pr checks {PR#} --watch --fail-fast
```

### Rules

1. Every required check green. No exceptions.
2. On failure: dispatch triage (see **Triage routing**), apply the `fix-hint`,
   commit, push, watch again. Expect shared triage output:
   `category`, `file`, `cause`, `fix-hint`, `handoff`.
3. Call something "flaky" only after two failures on the **same SHA**
   (`gh run rerun {run-id} --failed` once). One fail then one pass → proceed with
   a PR-body note. Two fails → not flaky; keep fixing. Prefer `test-triage` for
   that call — do not declare flake from the main session alone.
4. **Never** `--admin` or other branch-protection bypasses.

### Triage routing

Pick the agent from the failing check/log. When unclear, start with `log-triage`
and honor its `handoff`.

| Failure looks like… | Dispatch |
|---------------------|----------|
| Raw CI/job log, crash, setup, infra | `.cursor/agents/log-triage.md` |
| Failed assertion / test runner (pytest, jest, go test, …) | `.cursor/agents/test-triage.md` |
| Typechecker (tsc, mypy, pyright, …) | `.cursor/agents/type-triage.md` |
| Unclear | `log-triage` first → follow `handoff`: `test-triage` \| `type-triage` \| `none` |

Do not deep-parse stack traces in the main orchestrator session when a triage
agent can own them.

Also confirm local/project `verify` expectations from the `verify` skill if CI
does not cover them — do not claim green if a required gate was skipped.

---

## Phase 6: Specs check

```bash
git diff --name-only origin/{base}...HEAD
```

If the branch has **code** changes under the product tree and **no** `.specs/`
updates, stop and require a `docs:` specs update (or an explicit user waiver
recorded in the PR body). Pure docs/chore PRs may skip with a note.

---

## Phase 7: Merge

Only when **all** are true:

- [ ] Wave 1/2 verdicts are `APPROVE` or `APPROVE_WITH_CHANGES` (changes done)
- [ ] No unresolved CRITICAL or HIGH findings
- [ ] CI checks green
- [ ] Specs check satisfied or waived in the PR body
- [ ] Pre-flight mergeability still clean

```bash
gh pr view {PR#} --json mergeable,mergeStateStatus
# CONFLICTING / DIRTY → stop; do not merge
gh pr ready {PR#}
gh pr merge {PR#} --squash --auto
```

**Squash subject:** prefer Conventional Commits, ≤ 100 chars. Omitting `--body`
and leaving detail in the PR description is fine.

Confirm:

```bash
gh pr view {PR#} --json state,mergedAt,mergeCommit
```

---

## Phase 8: Cleanup

From the primary worktree (not necessarily the issue worktree):

```bash
git fetch origin
git checkout {base}
git pull
# If this work lived in a worktree:
git worktree remove .trees/{dir} --force   # when applicable
git branch -d {branch-name}                # after remote branch is gone / merged
```

Skip worktree removal if the user is still inside that worktree — tell them to
leave it first.

---

## Dry-run mode

`merge-pr --dry-run` (or user asks for review only):

- Run Phases 1–4 (and 3.5 if useful)
- Post/report findings
- Skip CI watch-to-merge, ready, merge, and cleanup

---

## Scope discipline

This skill merges PRs. It does **not**:

- Refactor beyond addressing findings
- Expand PR scope with unrelated changes
- Turn large review suggestions into merge blockers (file follow-ups instead)

---

## Handoff after merge

A merged PR is a hard task boundary. Prefer a fresh session.

1. Read the PR body for a `Milestone:` path under `.specs/milestones/` if present.
2. If a milestone remains with open issues → suggest continuing that milestone via
   `build` (or `plan-ready` if the next slice is not planned).
3. If unclear → ask; do not invent the next project.

### Output template

```markdown
**Merged**: PR #{PR#} → `{base}` (squash).
**Crew**: {summary of verdicts}.
**Follow-ups**: {list or none}.
**Cleanup**: {worktree/branch status}.

**Resume** (optional, after clearing context):
> Continue milestone `{slug}` from `.specs/milestones/{slug}/` — next open issue: …
```
