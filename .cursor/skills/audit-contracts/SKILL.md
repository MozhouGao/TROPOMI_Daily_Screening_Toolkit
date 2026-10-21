---
name: audit-contracts
description: >-
  Advisory schema/API contract audit. Use when reviewing public interfaces,
  validators, or compatibility. Does not merge or replace merge-pr.
disable-model-invocation: true
---

# audit-contracts

Advisory review of public contracts: APIs, schemas, validators, serialized
formats, and compatibility expectations.

## Phases

1. Fix scope (diff / paths / named contracts).
2. Inventory public interfaces touched or nominated.
3. Check against the lens below.
4. Emit findings (do not auto-fix unless asked).

## Lens

- Contracts are explicit (schema, types, OpenAPI, proto, …) — not only tribal knowledge
- Validation sits at trust / IO boundaries
- Compatibility: breaking changes acknowledged (versioning, migration, callers)
- Tests or fixtures lock the contract at the boundary
- Error shapes and empty/null cases are defined where consumers rely on them
- Docs/specs pointers exist for normative contracts (cite, don’t restate)

## Out of scope

| Concern | Where |
|---------|--------|
| Internal layering without a published contract | `audit-architecture` |
| General unit-test adequacy away from boundaries | `audit-tests` |
| Secret handling | security lane |

## Output

```text
id:       audit-contracts:{n}
severity: CRITICAL | HIGH | MEDIUM | LOW
file:     path (contract or producer/consumer)
summary:  …
rule:     optional (breaking, validation, versioning, follow-up, …)
fix:     direction only unless user asks to implement
```

Call out any silent breaking change as CRITICAL/HIGH.
