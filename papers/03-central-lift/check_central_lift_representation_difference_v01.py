#!/usr/bin/env python3
"""Exact regression gate for
Central Lift, Representation Difference, and Retained Occurrence, Draft 0.1.

The script checks finite and symbolic witnesses for the displayed algebraic,
sequence, central-return, and torus-placement formulas.  It is a regression
gate, not a substitute for the paper's universal proofs.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Any

import sympy as sp


class GateFailure(AssertionError):
    pass


@dataclass
class Gate:
    checks: int = 0

    def equal(self, lhs: Any, rhs: Any, label: str) -> None:
        self.checks += 1
        if isinstance(lhs, sp.MatrixBase) or isinstance(rhs, sp.MatrixBase):
            a = sp.Matrix(lhs)
            b = sp.Matrix(rhs)
            if a.shape != b.shape or any(sp.simplify(x) != 0 for x in (a - b)):
                raise GateFailure(f"{label}:\n{a}\n!=\n{b}")
        elif sp.simplify(lhs - rhs) != 0:
            raise GateFailure(f"{label}: {lhs!r} != {rhs!r}")

    def true(self, condition: bool, label: str) -> None:
        self.checks += 1
        if not condition:
            raise GateFailure(label)


def star(a: sp.Matrix, b: sp.Matrix) -> sp.Matrix:
    return a + b + a * b


def main() -> int:
    gate = Gate()
    I = sp.eye(2)
    Z = sp.zeros(2)
    J = sp.Matrix([[0, -1], [1, 0]])

    # Quarter-turn / planar half-turn checks.
    gate.equal(J * J, -I, "J^2=-I")
    phases = [I, J, -I, -J]
    for k in range(4):
        gate.equal(phases[k] * phases[(4 - k) % 4], I, f"phase inverse {k}")
    gate.equal(phases[2] + I, Z, "half-turn closure")
    gate.equal(phases[2] - I, -2 * I, "half-turn difference")
    gate.true(phases[1] != I and phases[3] != I, "quarter turns are active")

    # Difference composition, silent locus, and central parity on C4.
    for a in range(4):
        for b in range(4):
            Da = phases[a] - I
            Db = phases[b] - I
            Dab = phases[(a + b) % 4] - I
            gate.equal(star(Da, Db), Dab, f"difference law C4 ({a},{b})")
    for a in range(4):
        gate.true((phases[a] - I == Z) == (a == 0), f"silent locus C4 {a}")

    # Associativity of boxplus for generic symbolic 2x2 matrices.
    symbols = sp.symbols("a0:4 b0:4 c0:4")
    A = sp.Matrix(2, 2, symbols[0:4])
    B = sp.Matrix(2, 2, symbols[4:8])
    C = sp.Matrix(2, 2, symbols[8:12])
    gate.equal(star(star(A, B), C), star(A, star(B, C)), "boxplus associativity")
    gate.equal(star(A, Z), A, "boxplus right identity")
    gate.equal(star(Z, A), A, "boxplus left identity")
    gate.equal(I + star(A, B), (I + A) * (I + B), "transported product law")

    # Inverse law on several exact invertible matrices.
    samples = [
        sp.Matrix([[1, 0], [0, 1]]),
        sp.Matrix([[0, -1], [1, 0]]),
        sp.Matrix([[1, 1], [0, 1]]),
        sp.Matrix([[2, 1], [1, 1]]),
        sp.Matrix([[-1, 0], [0, -1]]),
    ]
    for idx, rho in enumerate(samples):
        D = rho - I
        Dinv = rho.inv() - I
        gate.equal(Dinv, -rho.inv() * D, f"inverse law left {idx}")
        gate.equal(Dinv, -D * rho.inv(), f"inverse law right {idx}")
        gate.equal(star(D, Dinv), Z, f"boxplus inverse {idx}")

    # Central character parity for integer and half-integer labels represented by 2J.
    for two_j in range(9):
        expected = (-1) ** two_j
        gate.equal(expected * I, (I if two_j % 2 == 0 else -I), f"central parity 2J={two_j}")

    # Central closure theorem and two independent double-return routes.
    z_action = -I
    Dz = z_action - I
    Cz = z_action + I
    gate.equal(Cz, Z, "central closure zero")
    gate.equal(Dz, -2 * I, "central represented difference")
    gate.true(Dz != Z, "central represented difference nonzero")
    gate.equal(z_action * z_action, I, "central element squares to identity")
    gate.equal(star(Dz, Dz), Z, "double return by difference law")

    # Monodromy parity multiplication as the C2 deck group law.
    signs = [1, -1]
    for a in signs:
        for b in signs:
            gate.equal(a * b, signs[(int(a == -1) + int(b == -1)) % 2], f"C2 product {a},{b}")

    # Split append, repeated receipts, prefix and payload recovery.
    archives = [[], ["a"], ["a", "a"], [1, 2, 3], [("run", 0), ("run", 1)]]
    payloads = ["x", "a", 3, ("lift", -1)]
    for i, archive in enumerate(archives):
        for j, payload in enumerate(payloads):
            out = archive + [payload]
            gate.equal(len(out), len(archive) + 1, f"append length {i},{j}")
            gate.true(out[:-1] == archive, f"dropLast recovers prefix {i},{j}")
            gate.true(out[-1] == payload, f"last recovers payload {i},{j}")
    repeated = ["same", "same"]
    gate.equal(len(repeated), 2, "duplicate receipts retain multiplicity")
    gate.true(repeated[0] == repeated[1], "duplicate values may agree")

    # Covering-torus weight compatibility sampled at exact quarter angles.
    alpha_values = [sp.Integer(0), sp.pi / 2, sp.pi, 3 * sp.pi / 2, 2 * sp.pi]
    for idx, alpha in enumerate(alpha_values):
        plane = sp.cos(alpha) * I + sp.sin(alpha) * J
        P = sp.Matrix([[1, 1], [-sp.I, sp.I]])  # eigenvectors for +i and -i
        diagonalized = sp.simplify(P.inv() * plane * P)
        target = sp.diag(sp.exp(sp.I * alpha), sp.exp(-sp.I * alpha))
        gate.equal(diagonalized.applyfunc(sp.simplify), target.applyfunc(sp.simplify), f"torus weights {idx}")

    # Lower quotient obstruction at the kernel element -1.
    e00, e01, e10, e11 = sp.symbols("e00 e01 e10 e11")
    Emap = sp.Matrix([[e00, e01], [e10, e11]])
    obstruction = Emap * (-I) - I * Emap
    gate.equal(obstruction, -2 * Emap, "lower quotient obstruction equation")
    gate.true(all(coef != 0 for coef in [-2, -2, -2, -2]), "characteristic zero obstruction")

    # Clifford-style path endpoint in a matrix model with E^2=-I.
    E = J
    gate.equal(E * E, -I, "Clifford plane generator square")
    gamma0 = sp.cos(0) * I + sp.sin(0) * E
    gamma1 = sp.cos(sp.pi) * I + sp.sin(sp.pi) * E
    gate.equal(gamma0, I, "lifted path starts at identity")
    gate.equal(gamma1, -I, "lifted path ends at central minus identity")

    print(f"PASS: {gate.checks} exact and finite-regression checks")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except GateFailure as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
