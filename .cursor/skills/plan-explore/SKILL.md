---
name: plan-explore
description: >-
  Explore planning options with roundtable personas before locking a plan.
  Use when approaches are ambiguous or the user asks to explore / pressure-test a plan.
disable-model-invocation: true
---

# plan-explore

Fan out roundtable personas, synthesize, update the milestone.

## Phases

1. Load the current plan (or draft one lightly under `.specs/milestones/{slug}/`).
2. Dispatch **all five** personas in parallel when possible:
   - `.cursor/agents/roundtable-minimalist.md`
   - `.cursor/agents/roundtable-skeptic.md`
   - `.cursor/agents/roundtable-pragmatist.md`
   - `.cursor/agents/roundtable-ops.md`
   - `.cursor/agents/roundtable-methodology.md`
3. Expect each to return the shared shape:
   `persona`, `stance` (`support` \| `concern` \| `block`), `points`, `asks`.
4. **Synthesize** (orchestrator owns this — personas do not merge themselves):
   - Agreements
   - Tensions (conflicting `points` / stances)
   - Union of `asks` that must be answered before build
   - Recommended path (and what to defer)
5. Update `.specs/milestones/{slug}/` with an exploration notes section (or linked doc).
6. Continue to `plan-ready`, or revise via `plan` / `decide` if `block` stances or
   critical asks remain.

## Synthesis rules

- Any `stance: block` → do not pretend the plan is ready; surface it for `plan` / `decide` / `plan-ready` **NO-GO**.
- Prefer cuts from `minimalist` when tensions are “size vs speed.”
- Record unresolved `asks` explicitly in the milestone — silent gaps are failures.
