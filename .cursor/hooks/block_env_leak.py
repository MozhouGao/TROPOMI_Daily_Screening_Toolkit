#!/usr/bin/env python3
"""Warn/deny commands that look like they dump secrets. Skeleton."""
from __future__ import annotations

import json
import re
import sys

LEAK = re.compile(
    r"(printenv|env\s*$|cat\s+[^\n]*\.env|echo\s+\$\{?(API|SECRET|TOKEN|PASSWORD))",
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
    if LEAK.search(command):
        print(
            json.dumps(
                {
                    "permission": "ask",
                    "user_message": "This command may expose environment secrets.",
                    "agent_message": "Avoid dumping env/secrets into the transcript.",
                }
            )
        )
        return

    print(json.dumps({"permission": "allow"}))


if __name__ == "__main__":
    main()
