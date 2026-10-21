#!/usr/bin/env python3
"""Deny obviously destructive shell commands. Skeleton — tighten for your org."""
from __future__ import annotations

import json
import re
import sys

DENY = re.compile(
    r"(rm\s+-rf\s+/|git\s+push\s+.*--force|git\s+reset\s+--hard|mkfs\.|dd\s+if=)",
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
    if DENY.search(command):
        print(
            json.dumps(
                {
                    "permission": "deny",
                    "user_message": "Blocked potentially destructive command.",
                    "agent_message": "Hook blocked a destructive shell command. Use a safer alternative.",
                }
            )
        )
        return

    print(json.dumps({"permission": "allow"}))


if __name__ == "__main__":
    main()
