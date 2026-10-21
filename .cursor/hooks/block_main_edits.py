#!/usr/bin/env python3
"""Block file edits while on main/master. Skeleton."""
from __future__ import annotations

import json
import subprocess
import sys


def current_branch() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return out.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ""


def main() -> None:
    # Consume stdin so Cursor does not see a broken pipe
    sys.stdin.read()
    branch = current_branch()
    if branch in {"main", "master"}:
        print(
            json.dumps(
                {
                    "permission": "deny",
                    "user_message": f"Edits blocked on {branch}. Use a worktree/feature branch.",
                    "agent_message": "Constitution: never edit on main. Create/switch to an issue worktree first.",
                }
            )
        )
        return

    print(json.dumps({"permission": "allow"}))


if __name__ == "__main__":
    main()
