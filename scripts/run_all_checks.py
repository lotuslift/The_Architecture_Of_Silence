#!/usr/bin/env python3
"""Run every shipped regression checker from its own artifact directory."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKS = [
    ("Complete 2x2 operator algebra", ROOT / "papers/01-complete-2x2-operator-algebra", "check_complete_2x2_operator_algebra_v02.py"),
    ("Condition and Constraint", ROOT / "papers/02-condition-and-constraint", "check_condition_and_constraint_v01.py"),
    ("Central lift", ROOT / "papers/03-central-lift", "check_central_lift_representation_difference_v01.py"),
    ("Four-phase / eight-state relations", ROOT / "papers/04-four-phase-eight-state", "check_four_phase_eight_state_relations_v01.py"),
    ("Conserved transport and return", ROOT / "papers/05-conserved-transport-return", "check_conserved_transport_return_v01.py"),
]


def main() -> int:
    failures: list[str] = []
    for label, directory, script in CHECKS:
        print(f"\n=== {label} ===", flush=True)
        result = subprocess.run([sys.executable, script], cwd=directory, check=False)
        if result.returncode:
            failures.append(f"{label} (exit {result.returncode})")

    if failures:
        print("\nFAILED:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print("\nPASS: all five regression checkers completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
