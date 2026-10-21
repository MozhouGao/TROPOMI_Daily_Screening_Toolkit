---
name: verify
description: >-
  Run project verification gates (lint, typecheck, tests, build). Use when the
  user says verify, check gates, or before merge-pr. Fail closed if a required
  gate did not run.
disable-model-invocation: true
---

# verify

Run the project gate suite and report pass/fail per gate. Used by `build` and
before trusting a merge path; `merge-pr` still requires CI green separately.

## Configure (per product repo)

Replace placeholders with this repo’s real commands. Keep this table honest —
do not claim a gate you cannot run.

| Gate | Required | Command (edit me) |
|------|----------|-------------------|
| lint | yes | `pnpm --filter ws-emission-dashboard lint` |
| typecheck | yes | `pnpm --filter ws-emission-dashboard exec tsc --noEmit` |
| test | yes | `pnpm --filter ws-emission-dashboard test` |
| build | optional | `pnpm --filter ws-emission-dashboard build` |

Add/remove rows to match the product. Mark optional gates explicitly.

## Phases

1. Confirm cwd is the worktree/repo root you intend to verify.
2. Run each **required** gate (and optional ones when relevant to the change).
3. On failure: dispatch triage (see **Triage routing** below). Load persona
   prompts from `.cursor/agents/{name}.md` and run each agent with that file as
   its instructions plus the failing gate output.
4. Re-run failed gates after fixes until green or blocked.
5. Emit the report below.

### Triage routing

Pick the agent from the failure. When unclear, start with `log-triage` and honor
its `handoff`.

| Failure looks like… | Dispatch |
|---------------------|----------|
| Raw log, crash, setup, infra, or gate could not run | `.cursor/agents/log-triage.md` |
| Failed assertion / test runner (vitest, jest, pytest, …) | `.cursor/agents/test-triage.md` |
| Typechecker (`tsc`, mypy, pyright, …) | `.cursor/agents/type-triage.md` |
| Unclear | `log-triage` first → follow `handoff`: `test-triage` \| `type-triage` \| `none` |

Do not deep-parse stack traces in the main orchestrator session when a triage
agent can own them. Same routing table as `merge-pr` Phase 5.

## Output

```markdown
| Gate | Status | Notes |
|------|--------|-------|
| lint | PASS/FAIL/SKIP | |
| typecheck | PASS/FAIL/SKIP | |
| test | PASS/FAIL/SKIP | |
| build | PASS/FAIL/SKIP | |

**Overall:** GREEN | RED | BLOCKED
```

| Overall | Meaning |
|---------|---------|
| `GREEN` | All required gates PASS |
| `RED` | At least one required gate FAIL |
| `BLOCKED` | Could not run a required gate (missing tool, env, etc.) |

## Rules

1. **Fail closed:** do not claim `GREEN` if a required gate did not run.
2. `SKIP` only for optional gates or an explicit user waiver recorded in notes.
3. Do not use `--admin`, ignore scripts, or delete failing tests to get green.
4. Prefer synthetic/local fixtures; do not require production secrets for verify.

## Next

| Overall | Next |
|---------|------|
| `GREEN` | Continue `build` or proceed to `merge-pr` |
| `RED` | Fix via triage (`.cursor/agents/*-triage.md`) / `build`+`tdd` / code; re-verify |
| `BLOCKED` | Unblock environment or document waiver; do not merge as green |
