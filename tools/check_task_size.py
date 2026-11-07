#!/usr/bin/env python3
"""Check the aggregate task budget against the task-base Git reference."""

import argparse
from pathlib import Path
import subprocess
import sys

LIMIT = 10_000_000


def git(*args):
    return subprocess.check_output(["git", *args])


def check(base="task-base", staged=False):
    root = Path(git("rev-parse", "--show-toplevel").decode().strip())
    import os
    os.chdir(root)
    git("rev-parse", "--verify", base + "^{commit}")
    options = ["--cached"] if staged else []
    patch = git("diff", "--binary", "--no-ext-diff", "--no-textconv", *options, base, "--")
    changed = git("diff", "--name-only", "--no-renames", "-z", *options, base, "--").split(b"\0")
    tracked_size = 0
    for raw_path in filter(None, changed):
        path = os.fsdecode(raw_path)
        if staged:
            entry = git("ls-files", "--stage", "--", path).decode().strip()
            if not entry:
                continue  # Deletion
            mode, oid, stage = entry.split("\t", 1)[0].split()
            if mode == "160000":
                raise ValueError("Nested repositories/submodules must be submitted separately.")
            tracked_size += int(git("cat-file", "-s", oid))
        else:
            file = root / path
            if file.is_symlink():
                tracked_size += len(os.fsencode(os.readlink(file)))
            elif file.is_file():
                tracked_size += file.stat().st_size
            elif file.is_dir():
                raise ValueError("Nested repositories/submodules must be submitted separately.")
    untracked_size = 0
    if not staged:
        for raw_path in filter(None, git("ls-files", "--others", "--exclude-standard", "-z").split(b"\0")):
            file = root / os.fsdecode(raw_path)
            if file.is_symlink():
                untracked_size += len(os.fsencode(os.readlink(file)))
            elif file.is_file():
                untracked_size += file.stat().st_size
            elif file.is_dir():
                raise ValueError("Untracked nested repositories must be submitted separately.")
    # Include both the binary patch and changed file payloads as a conservative budget.
    total = len(patch) + tracked_size + untracked_size
    print(f"Task budget: {total:,} / {LIMIT:,} bytes")
    print(f"Binary patch: {len(patch):,}; changed files: {tracked_size:,}; untracked files: {untracked_size:,}")
    if total >= LIMIT:
        print("FAIL: task changes must stay below 10 MB (decimal).", file=sys.stderr)
        return 1
    print("PASS: below 10 MB. Ignored local artifacts are not part of the submission.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="task-base")
    parser.add_argument("--staged", action="store_true", help="Check the index, as used by the commit hook.")
    args = parser.parse_args()
    try:
        return check(args.base, args.staged)
    except (subprocess.CalledProcessError, ValueError, OSError) as exc:
        print(f"Cannot verify the task budget: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
