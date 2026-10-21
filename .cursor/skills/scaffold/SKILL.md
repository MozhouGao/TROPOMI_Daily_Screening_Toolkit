---
name: scaffold
description: >-
  Create file/package skeletons from an approved plan without implementing behavior.
  Use when the user says scaffold, stub, or bootstrap structure.
disable-model-invocation: true
---

# scaffold

From an approved (or clearly scoped) milestone plan, create empty/minimal
structure only — no business logic.

## Preconditions

- Prefer a **GO** from `plan-ready`, or an explicit user request to scaffold early
- Work on a worktree/branch — not `main`
- Know the target paths from `.specs/milestones/{slug}/` or the user

## Phases

1. Read the plan’s approach + issues; list files/packages to create.
2. Match **existing** repo conventions (layout, naming, exports, test placement).
3. Create stubs only: empty modules, interfaces, placeholders, `TODO(build):`
   markers tied to acceptance criteria where useful.
4. Optionally open/update a **draft** PR if the user wants the scaffold landed.
5. Stop. Hand off to `build` (with `tdd` for behavior).

## Rules

- No business logic, no “temporary” real implementations.
- No drive-by refactors of unrelated code.
- Do not invent new frameworks or folder schemes unless the plan/ADR says so.
- Keep changes reviewable and draft-PR friendly.
- If structure choices are still undecided → `decide` or `plan`, not guesswork.

## Output contract

- [ ] Listed stub paths created
- [ ] Conventions matched
- [ ] TODOs point at plan criteria / issues where helpful
- [ ] No behavior claimed as done

## Next

`build` (TDD slices for testable criteria) → `verify` → `merge-pr`.
