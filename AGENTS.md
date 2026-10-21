# Agent Framework Constitution

Framework-only SDLC harness for Cursor. Replace bracketed placeholders when you adopt this in a product repo.

## Constitutional rules

### 1. Always work in a worktree

Never implement on `main`. All planning and coding happens on an issue branch / worktree.

```bash
git fetch origin
git worktree add .trees/{ISSUE-ID} -b issue/{ISSUE-ID}-{slug}
cd .trees/{ISSUE-ID}
```

### 2. Always commit and push

Meaningful changes are committed and pushed. Uncommitted work is invisible; unpushed work is at risk.

### 3. Always use draft pull requests

Never push directly to `main`. Create PRs as **draft**. Never use `--admin` to bypass branch protection.

### 4. Always update specs before merge

Discoveries, decisions, and edge cases flow into `.specs/` on the same branch (`docs:` commits).

### 5. Canonical review at merge-pr

`build` leaves PRs as drafts. `merge-pr` is the only review-crew dispatch point.

Wave 1 (parallel): `reviewer`, `architect`, `security` (when relevant)  
Wave 2 (last): `adversary`

Each crew agent returns a verdict (`APPROVE` | `APPROVE_WITH_CHANGES` | `BLOCK` | `SKIP`) and id’d findings (`{agent}:{n}`). Address CRITICAL and HIGH before marking the PR ready.

### 6. Always write tests for new code

Tests are part of done. Prefer deterministic, synthetic fixtures over production snapshots. For new testable behavior in `build`, prefer `tdd` slices.

### 7. Prefer simplicity

Propose the smallest change that satisfies the requirement. No premature platforms.

### 8. Never commit secrets

Credentials, API keys, tokens, and customer-identifying data must not enter git or session transcripts.

---

## Skill catalog

| Skill | When to use | Dispatches |
|-------|-------------|------------|
| `plan` | Scope a milestone into `.specs/milestones/` | — |
| `decide` | Capture an ADR in `.specs/decisions/` | — |
| `plan-explore` | Pressure-test / explore options before locking a plan | All five `roundtable-*` |
| `plan-ready` | Gate: is the plan ready to build? | Roundtable when gaps remain or user asks (else reuse explore notes) |
| `scaffold` | File/package skeleton from a plan (no behavior) | — |
| `build` | Implement on a **draft** PR until gates are green | `tdd` per testable acceptance criterion |
| `verify` | Run project verification gates (fail closed) | Triage on fail (`log-triage` / `test-triage` / `type-triage`) |
| `merge-pr` | Review crew → CI triage → ready → squash-merge | Wave 1–2 crew; triage on CI fail |
| `audit` | Advisory audit orchestrator | Routes to audit-* skills |
| `audit-architecture` | Architecture / boundary audit | — |
| `audit-contracts` | Schema / API contract audit | — |
| `audit-tests` | Test adequacy audit | — |

Skills live under `.cursor/skills/*/SKILL.md`. A Kimi Code CLI mirror lives under `.kimi-code/` (regenerate with `scripts/sync-kimi-from-cursor.py`).

## Agent catalog

Persona prompts live in `.cursor/agents/`. Skills invoke them; Cursor does not auto-load that folder.

### Review crew (`merge-pr`)

| Agent | Role | Output |
|-------|------|--------|
| `reviewer` | Code quality on the diff | Verdict + `reviewer:{n}` findings |
| `architect` | Placement, layering, coupling | Verdict + `architect:{n}` findings |
| `security` | Secrets / trust boundaries (skip if not triggered) | Verdict + `security:{n}` findings, or `SKIP` |
| `adversary` | Wave 2 meta-review; requires Wave 1 ids | Verdict + `adversary:{n}` + `crossChecks` |

### CI triage (`merge-pr` Phase 5)

| Agent | When | Output |
|-------|------|--------|
| `log-triage` | Unclear CI / crash / setup / infra (default entry) | `category`, `file`, `cause`, `fix-hint`, `handoff` |
| `test-triage` | Failed assertions / test runner | same |
| `type-triage` | Typechecker failures | same |

Routing table: see `merge-pr` → **Triage routing**.

### Build helper (`build`)

| Agent | When | Output |
|-------|------|--------|
| `tdd` | Each testable acceptance criterion (new/changed behavior) | `criterion`, `test`, `status` (`red` \| `green` \| `blocked`), `notes` |

### Roundtable (`plan-explore`, optionally `plan-ready`)

| Agent | Lane | Output |
|-------|------|--------|
| `roundtable-minimalist` | Cut / defer scope | `persona`, `stance`, `points`, `asks` |
| `roundtable-skeptic` | Assumptions / failure modes | same |
| `roundtable-pragmatist` | Sequencing / speed vs risk | same |
| `roundtable-ops` | Rollout / observe / rollback / support | same |
| `roundtable-methodology` | Metrics / validation / DoD | same |

`stance`: `support` | `concern` | `block`. Any `block` or critical open `asks` → plan is not ready.

## Work loops

**Feature / fix**

```
plan → (plan-explore) → plan-ready → (scaffold) → build[+tdd] → verify → merge-pr
```

**Decision only**

```
decide → (optional) plan-ready
```

**Advisory**

```
audit | audit-architecture | audit-contracts | audit-tests
```
