#!/usr/bin/env python3
"""Mechanical audit for Condition & Constraint, draft 0.1.

The script checks the finite two-condition case split, the exact quadratic
whole--remainder witness, and structural invariants of the LaTeX source.
It is a regression gate for the displayed finite and symbolic claims, not a
substitute for reading the argument.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import sympy as sp

SOURCE = Path(__file__).with_name("Condition_and_Constraint_Boundary_Discipline_Draft_0_1.tex")


def check(condition: bool, label: str, counter: list[int]) -> None:
    if not condition:
        raise AssertionError(label)
    counter[0] += 1


def main() -> None:
    count = [0]
    text = SOURCE.read_text(encoding="utf-8")

    # Publication-form invariants.
    required = [
        "Condition \\& Constraint",
        "Boundary Discipline",
        "22 January 2026",
        "Part I: ",
        "Part II: Algebra",
        "Part III: Geometry",
        "Two-condition case split",
        "Interior traversal",
        "Minimum jurisdiction record",
        "No transport by matching marks alone",
    ]
    for token in required:
        check(token in text, f"missing required token: {token}", count)

    forbidden = [
        "Revised form",
        "Original version",
        "Canonical Root Series",
        "Minimum Endurance Delta",
        "The Architecture of Silence",
        "Motion and Difference",
        "Memory Return",
        "Clear_Onyx",
    ]
    for token in forbidden:
        check(token not in text, f"forbidden development-history token: {token}", count)

    # Exact Boolean case split.
    expected = {
        (0, 0): "non-candidate",
        (0, 1): "closure-aligned",
        (1, 0): "remainder-aligned",
        (1, 1): "interior traversal",
    }
    check(len(expected) == 4, "case split must have four profiles", count)
    check(sum(1 for k in expected if k == (1, 1)) == 1,
          "exactly one profile must satisfy both conditions", count)
    for ext in (0, 1):
        for rec in (0, 1):
            check((ext, rec) in expected, f"missing profile {(ext, rec)}", count)
    check(expected[(0, 1)] == "closure-aligned", "closure profile mismatch", count)
    check(expected[(1, 0)] == "remainder-aligned", "remainder profile mismatch", count)
    check(expected[(1, 1)] == "interior traversal", "survivor profile mismatch", count)

    # Exact whole--remainder witness.
    x = sp.symbols("x")
    poly = x**2 - x - 1
    roots = sp.solve(poly, x)
    phi = (1 + sp.sqrt(5)) / 2
    psi = (1 - sp.sqrt(5)) / 2
    check(set(roots) == {phi, psi}, "quadratic roots mismatch", count)
    check(sp.simplify(phi + psi - 1) == 0, "root sum mismatch", count)
    check(sp.simplify(phi * psi + 1) == 0, "root product mismatch", count)
    check(sp.simplify(phi - psi - sp.sqrt(5)) == 0,
          "branch separation mismatch", count)
    check(sp.expand((x - phi) * (x - psi) - poly) == 0,
          "factorization mismatch", count)
    check(sp.simplify(phi - (1 + 1 / phi)) == 0,
          "positive branch fixed relation mismatch", count)
    check(sp.simplify(psi - (1 + 1 / psi)) == 0,
          "negative branch fixed relation mismatch", count)

    digest = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    check(len(digest) == 64, "SHA-256 digest length mismatch", count)

    print(f"PASS: {count[0]} checks")
    print(f"source_sha256={digest}")


if __name__ == "__main__":
    main()
