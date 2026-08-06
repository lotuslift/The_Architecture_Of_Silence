#!/usr/bin/env python3
"""Exact regression gate for Four-Phase Addressing and Eight-State Relations.

The script checks finite instances of the paper's displayed identities. It is a
regression gate, not a substitute for the universal proofs in the paper.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from typing import Dict, Iterable, List, Sequence, Tuple

CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Four-phase / occurrence addressing
# ---------------------------------------------------------------------------

def time_read(n: int) -> Tuple[int, int]:
    return n % 4, n // 4


def unread(p: int, c: int) -> int:
    return 4 * c + p


for n in range(2000):
    p, c = time_read(n)
    check(0 <= p < 4, f"phase out of range at n={n}")
    check(unread(p, c) == n, f"unread(time_read({n})) failed")
    p2, c2 = time_read(n + 1)
    expected = (p + 1, c) if p < 3 else (0, c + 1)
    check((p2, c2) == expected, f"successor carry failed at n={n}")
    check(time_read(n + 4) == (p, c + 1), f"four-step return failed at n={n}")

for p in range(4):
    for c in range(500):
        check(time_read(unread(p, c)) == (p, c), f"time_read(unread({p},{c})) failed")

run_id = "gate-run"
for n in range(1000):
    occurrence = (run_id, n)
    occurrence_4 = (run_id, n + 4)
    check(occurrence != occurrence_4, f"occurrence erased at n={n}")

# Two-window address
addresses = {(ell, p): 4 * ell + p for ell in (0, 1) for p in range(4)}
check(set(addresses.values()) == set(range(8)), "two-window address is not bijective")
check(len(set(addresses.values())) == 8, "two-window addresses are not distinct")

# ---------------------------------------------------------------------------
# Pair count and alternating relation record
# ---------------------------------------------------------------------------
vertices = tuple(range(8))
unordered_edges = tuple(combinations(vertices, 2))
check(len(unordered_edges) == 28, "K8 must have 28 unordered edges")
check(2 * len(unordered_edges) == 56, "K8 must have 56 oriented edges")

# Deterministic skew relation record.
d: Dict[Tuple[int, int], int] = {}
for i, j in unordered_edges:
    value = 3 * i - 2 * j + 5
    d[(i, j)] = value
    d[(j, i)] = -value
for i in vertices:
    d[(i, i)] = 0

for i in vertices:
    check(d[(i, i)] == 0, f"diagonal relation not zero at {i}")
    for j in vertices:
        check(d[(j, i)] == -d[(i, j)], f"antisymmetry failed at {(i, j)}")

# Additive triangle criterion: exhaustive forward direction on 3^8 potentials.
for s in product((-1, 0, 1), repeat=8):
    for i, j, k in combinations(vertices, 3):
        dij = s[j] - s[i]
        djk = s[k] - s[j]
        dki = s[i] - s[k]
        check(dij + djk + dki == 0, f"potential triangle did not close for {s}")
    reconstructed = tuple((s[i] - s[0]) for i in vertices)
    check(reconstructed[0] == 0, "reconstruction base must be zero")
    for i, j in unordered_edges:
        check(reconstructed[j] - reconstructed[i] == s[j] - s[i], "potential reconstruction failed")

# Exhaustive converse on K4 with edge values in {-1,0,1}.
k4_vertices = range(4)
k4_edges = list(combinations(k4_vertices, 2))
for values in product((-1, 0, 1), repeat=len(k4_edges)):
    dd: Dict[Tuple[int, int], int] = {}
    for edge, value in zip(k4_edges, values):
        i, j = edge
        dd[(i, j)] = value
        dd[(j, i)] = -value
    for i in k4_vertices:
        dd[(i, i)] = 0
    closed = all(
        dd[(i, j)] + dd[(j, k)] + dd[(k, i)] == 0
        for i, j, k in combinations(k4_vertices, 3)
    )
    if closed:
        s = [dd[(0, i)] for i in k4_vertices]
        for i, j in k4_edges:
            check(dd[(i, j)] == s[j] - s[i], "closed K4 relation failed to reconstruct")

# ---------------------------------------------------------------------------
# K8 one-factorization and scanner
# ---------------------------------------------------------------------------
INF = "inf"
Z7 = tuple(range(7))


def unordered(a: object, b: object) -> frozenset:
    return frozenset((a, b))


def matching(a: int) -> Tuple[frozenset, ...]:
    return (
        unordered(INF, a),
        unordered((a + 1) % 7, (a - 1) % 7),
        unordered((a + 2) % 7, (a - 2) % 7),
        unordered((a + 3) % 7, (a - 3) % 7),
    )

matchings = {a: matching(a) for a in Z7}
all_vertices = {INF, *Z7}
seen_edges: List[frozenset] = []
for a, edges in matchings.items():
    check(len(edges) == 4, f"matching {a} does not have four edges")
    covered = set().union(*edges)
    check(covered == all_vertices, f"matching {a} is not perfect")
    check(sum(len(edge) for edge in edges) == 8, f"matching {a} has repeated vertex")
    seen_edges.extend(edges)

check(len(seen_edges) == 28, "factorization does not list 28 edges")
check(len(set(seen_edges)) == 28, "factorization repeats an edge")
complete_k8 = {
    unordered(x, y)
    for x, y in combinations((INF, *Z7), 2)
}
check(set(seen_edges) == complete_k8, "factorization does not cover K8")


def kappa(a: int, b: int) -> frozenset:
    if b == 0:
        return unordered(INF, a)
    return unordered((a + b) % 7, (a - b) % 7)


def scan(q: int) -> frozenset:
    return kappa(q % 7, q % 4)

scanner_edges = [scan(q) for q in range(28)]
check(len(set(scanner_edges)) == 28, "scanner is not injective")
check(set(scanner_edges) == complete_k8, "scanner is not surjective")
for q in range(28):
    a, b = q % 7, q % 4
    check(scan(q) in matchings[a], f"scanner round mismatch at q={q}")
    check(scan(q) == kappa(a, b), f"scanner coordinate mismatch at q={q}")

# ---------------------------------------------------------------------------
# Relative phase over exact finite roots of unity
# ---------------------------------------------------------------------------
MOD = 4  # exact exponent model of the subgroup mu_4 <= U(1)


def pair_read(exponents: Sequence[int], i: int, j: int) -> int:
    """Exponent of z_j * conjugate(z_i) modulo MOD."""
    return (exponents[j] - exponents[i]) % MOD

relative_representatives = set()
pair_records = set()
for z in product(range(MOD), repeat=8):
    canonical = tuple((z[j] - z[0]) % MOD for j in range(1, 8))
    relative_representatives.add(canonical)

    # Invariance under all common phase shifts.
    base_record = tuple(pair_read(z, i, j) for i, j in unordered_edges)
    pair_records.add(base_record)
    for lam in range(MOD):
        shifted = tuple((x + lam) % MOD for x in z)
        shifted_record = tuple(pair_read(shifted, i, j) for i, j in unordered_edges)
        check(shifted_record == base_record, "common-phase invariance failed")

    # Pair identities.
    for i in vertices:
        check(pair_read(z, i, i) == 0, "multiplicative diagonal identity failed")
        for j in vertices:
            check((pair_read(z, i, j) + pair_read(z, j, i)) % MOD == 0, "inverse relation failed")
    for i, j, k in combinations(vertices, 3):
        check(
            (pair_read(z, i, j) + pair_read(z, j, k) - pair_read(z, i, k)) % MOD == 0,
            "multiplicative composition failed",
        )
        check(
            (pair_read(z, i, j) + pair_read(z, j, k) + pair_read(z, k, i)) % MOD == 0,
            "multiplicative triangle failed",
        )

check(len(relative_representatives) == MOD ** 7, "quotient must have MOD^7 relative states")
check(len(pair_records) == MOD ** 7, "pair record must determine relative state")

# Converse reconstruction for every exact relative state in mu_4^7.
for rel in product(range(MOD), repeat=7):
    z = (0, *rel)
    record = {(i, j): pair_read(z, i, j) for i in vertices for j in vertices}
    reconstructed = tuple(record[(0, j)] for j in range(1, 8))
    check(reconstructed == rel, "relative phase reconstruction failed")
    for i, j, k in combinations(vertices, 3):
        check((record[(i, j)] + record[(j, k)] + record[(k, i)]) % MOD == 0, "reconstructed triangle failed")

# ---------------------------------------------------------------------------
# Lightweight exact state-transition interface witness
# ---------------------------------------------------------------------------
Matrix = Tuple[Tuple[Fraction, Fraction], Tuple[Fraction, Fraction]]

I: Matrix = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))
Q: Matrix = ((Fraction(1), Fraction(1)), (Fraction(1), Fraction(0)))
C1: Matrix = ((Fraction(0), Fraction(-1)), (Fraction(1), Fraction(0)))
OMEGA: Matrix = ((Fraction(0), Fraction(1)), (Fraction(-1), Fraction(0)))
G: Matrix = ((Fraction(-2), Fraction(1)), (Fraction(1), Fraction(2)))


def mmul(a: Matrix, b: Matrix) -> Matrix:
    return (
        (
            a[0][0] * b[0][0] + a[0][1] * b[1][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1],
        ),
        (
            a[1][0] * b[0][0] + a[1][1] * b[1][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1],
        ),
    )


def madd(*terms: Matrix) -> Matrix:
    return tuple(
        tuple(sum(term[i][j] for term in terms) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def smul(c: Fraction, a: Matrix) -> Matrix:
    return tuple(tuple(c * a[i][j] for j in range(2)) for i in range(2))  # type: ignore[return-value]


def nf(a: Matrix) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    aa, b = a[0]
    c, dd = a[1]
    return (
        Fraction(2 * aa - b - c + 3 * dd, 5),
        Fraction(aa + 2 * b + 2 * c - dd, 5),
        Fraction(b - c, 2),
        Fraction(-2 * aa + b + c + 2 * dd, 10),
    )


def evaluate(coords: Sequence[Fraction]) -> Matrix:
    return madd(smul(coords[0], I), smul(coords[1], Q), smul(coords[2], OMEGA), smul(coords[3], G))

state_matrix = I
state_phase = tuple(range(8))  # exponents mod 4, reduced below
ledger: List[dict] = []
word = ""
for n, generator in enumerate((Q, C1, Q, Q, C1, C1, Q, C1), start=1):
    symbol = "Q" if generator == Q else "C1"
    before = state_matrix
    state_matrix = mmul(state_matrix, generator)
    coords = nf(state_matrix)
    check(evaluate(coords) == state_matrix, f"normal-form interface failed at transition {n}")
    p, c = time_read(n)
    check(unread(p, c) == n, f"state address failed at transition {n}")
    # Equivariant exact channel update: common quarter-turn shift.
    state_phase = tuple((x + 1) % MOD for x in state_phase)
    before_pairs = tuple(pair_read(tuple((x - 1) % MOD for x in state_phase), i, j) for i, j in unordered_edges)
    after_pairs = tuple(pair_read(state_phase, i, j) for i, j in unordered_edges)
    check(after_pairs == before_pairs, f"common phase changed relative seams at transition {n}")
    word += symbol
    receipt = {
        "n": n,
        "generator": symbol,
        "word": word,
        "coords": coords,
        "phase": p,
        "cycle": c,
        "pair_count": len(after_pairs),
    }
    ledger.append(receipt)
    check(len(ledger) == n, f"ledger length failed at transition {n}")
    check(ledger[-1] == receipt, f"last receipt failed at transition {n}")
    check(before != state_matrix or symbol == "", f"operator did not update at transition {n}")

print(f"PASS: {CHECKS:,} exact and finite-regression checks")
