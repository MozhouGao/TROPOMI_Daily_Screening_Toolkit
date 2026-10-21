#!/usr/bin/env python3
"""Soft-check commit message style on `git commit`. Skeleton — Conventional Commits."""
from __future__ import annotations

import json
import re
import sys

COMMIT = re.compile(r"\bgit\s+commit\b", re.I)
# Detect -m "..." inline messages that are clearly not conventional
BAD_MSG = re.compile(
    r'-m\s+[\'"](?!(feat|fix|docs|style|refactor|perf|test|build|ci|chore|revert)(\(.+\))?!?:).+[\'"]',
    re.I,
)


def main() -> None:
    raw = sys.stdin.read()
    try:
        payload = json.loads(raw) if raw.strip() else {}
    except json.JSONDecodeError:
        print(json.dumps({"permission": "allow"}))
        return

    command = payload.get("command") or ""
    if COMMIT.search(command) and BAD_MSG.search(command):
        print(
            json.dumps(
                {
                    "permission": "ask",
                    "user_message": "Commit message may not follow Conventional Commits.",
                    "agent_message": "Prefer: type(scope): subject  (feat|fix|docs|refactor|test|chore|…)",
                }
            )
        )
        return

    print(json.dumps({"permission": "allow"}))


if __name__ == "__main__":
    main()
