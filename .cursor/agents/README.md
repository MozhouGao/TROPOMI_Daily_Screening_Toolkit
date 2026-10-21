# Agent personas

Markdown prompts for review crew, triage, and roundtable personas.

Cursor does **not** auto-discover this folder the way Claude Code loads `.claude/agents/`. Skills should **read these files** and dispatch work with that persona’s instructions:

| Skill | Agents loaded from here |
|-------|-------------------------|
| `build` | `tdd.md` (per acceptance criterion) |
| `merge-pr` | `reviewer`, `architect`, `security`, `adversary` |
| `merge-pr` / `verify` (on fail) | `log-triage`, `test-triage`, `type-triage` |
| `plan-explore` / `plan-ready` | all five `roundtable-*.md` |

Customize checklists for your stack; keep names stable so skills can reference them.
