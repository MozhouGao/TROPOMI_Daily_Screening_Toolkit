---
name: plan-ready
description: >-
  Gate whether a plan is ready to build. Use before implementation when the user
  says plan-ready, ready to build, or asks if the plan is complete. May run or
  reuse roundtable pressure-test when gaps remain.
disable-model-invocation: true
---

# plan-ready

Readiness checklist → roundtable when needed → GO / NO-GO.

## Checklist

- [ ] Problem and scope are clear
- [ ] Acceptance criteria are testable
- [ ] Out of scope is explicit
- [ ] Risks / open questions listed
- [ ] Spec files exist under `.specs/`
- [ ] No unresolved blockers that require a decision (`decide` first)

## Roundtable

Run or reuse a pressure-test when **any** of:

- User asks to pressure-test / explore readiness
- Checklist items are weak or open questions remain
- Milestone has no prior exploration notes from `plan-explore`

**Reuse** existing `plan-explore` notes if they are fresh and cover the current
plan revision; otherwise dispatch all five `.cursor/agents/roundtable-*.md`
personas (same contract as `plan-explore`: `persona`, `stance`, `points`, `asks`).

### Roundtable → gate

| Signal | Effect |
|--------|--------|
| Any `stance: block` | **NO-GO** until addressed or explicitly waived by the user in the milestone |
| Unresolved critical `asks` | **NO-GO** (answer via `plan` / `decide`, then re-run) |
| Only `support` / minor `concern` | Checklist decides GO / NO-GO |

## Output

- **GO** — proceed to `scaffold` / `build`
- **NO-GO** — list gaps (checklist misses, `block` stances, open `asks`); return to `plan`, `decide`, or `plan-explore`
