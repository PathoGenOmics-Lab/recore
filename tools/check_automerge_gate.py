#!/usr/bin/env python3
"""Fail when the Dependabot auto-merge gate waits for a check that no job publishes.

The gate in .github/workflows/dependabot-auto-merge.yml names the checks it waits
for as `state_of '<workflow name>' '<check name>'`. A check name is what GitHub
publishes for a job: its `name:` if it has one, otherwise its id. When the two
drift apart nothing errors. The gate just never sees the check finish, waits out
its timeout on every pull request, and merges nothing. That happened here: it
waited for "test" while ci.yml published "Build, test, clippy, fmt".

Run from the repository root. Exits 1 and says which name is missing.
"""

import re
import sys
from pathlib import Path

import yaml

WORKFLOWS = Path(".github/workflows")
GATE = WORKFLOWS / "dependabot-auto-merge.yml"


def published_checks():
    """Map each workflow's name to the set of check names its jobs publish."""
    out = {}
    for path in sorted(WORKFLOWS.glob("*.yml")):
        doc = yaml.safe_load(path.read_text())
        name = doc.get("name", path.stem)
        jobs = doc.get("jobs") or {}
        out[name] = {str(job.get("name", job_id)) for job_id, job in jobs.items()}
    return out


def gate_expectations():
    text = GATE.read_text()
    pairs = re.findall(r"state_of '([^']+)' '([^']+)'", text)
    if not pairs:
        sys.exit(f"{GATE}: found no state_of '<workflow>' '<check>' calls to verify")
    return pairs


def main():
    published = published_checks()
    missing = []
    for workflow, check in gate_expectations():
        if check not in published.get(workflow, set()):
            have = sorted(published.get(workflow, set())) or "no such workflow"
            missing.append(f"  {workflow!r} / {check!r}: that workflow publishes {have}")
    if missing:
        print(f"{GATE} waits for checks that no job publishes:", file=sys.stderr)
        print("\n".join(missing), file=sys.stderr)
        sys.exit(1)
    print(f"{GATE}: every check the gate waits for is published")


if __name__ == "__main__":
    main()
