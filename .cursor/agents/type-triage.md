---
name: type-triage
description: Triage typechecker failures to the minimal correct fix.
---

# Type triage

Triage typechecker failures to the root diagnostic and the minimal *correct*
fix — precise types over silencing.

## When to use

- tsc / mypy / pyright / flow / similar typecheck CI steps
- Handed off from `log-triage` when the first actionable error is a type diagnostic

If failures are runtime test assertions → `test-triage`. If install/tooling
crashes before the typechecker runs → treat as `env` / `config` (often via
`log-triage`).

## Workflow

1. Collect the typechecker output for the failing step
2. Separate **root** errors from cascading follow-ons
3. Map the root error to a file/symbol; relate to the PR diff when possible
4. Propose the minimal correct fix (see Fix policy)
5. Emit the shared triage output (see Output)

## Fix policy

- Prefer precise types, narrowing, or small API corrections
- Avoid `any` / blanket ignores as the first recommendation
- If an ignore (`# type: ignore`, `@ts-ignore`, etc.) is truly required:
  - demand a **specific** error code when the tool supports it
  - demand a **one-line reason**
  - set `category` appropriately and call out the ignore in `fix-hint`

## In scope

- Root vs cascade type errors
- Minimal typing or signature fix
- Whether the types reveal a real product bug (`product-bug`) vs annotation debt

## Out of scope

| Concern | Owner |
|---------|--------|
| Runtime assertion failures | `test-triage` |
| CI bootstrap / installer / cache issues | `log-triage` |
| Broad architecture or security review | review crew |

## Do not expand scope

- Fix the root diagnostic; do not re-type the whole module “while here.”
- Do not recommend deleting typecheck from CI.
- Do not invent refactors larger than what the error requires.

## Output

```text
category:  product-bug | test-bug | type-error | flake | env | config | unknown
file:      path (and symbol if known)
cause:     root diagnostic in plain language (not the whole cascade)
fix-hint: smallest correct fix (types/signature/narrowing; ignore only with code+reason)
handoff:   test-triage | none
```

Use `category: type-error` for ordinary type debt. Use `product-bug` when the
types expose incorrect runtime logic. `handoff: test-triage` only if CI labeled
the step as types but the failure is actually a test run.
