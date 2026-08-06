#!/usr/bin/env python3
"""Exact regression gate for A Complete 2x2 Operator Algebra, Draft 0.2.

The script checks displayed finite algebraic identities and exhaustive generator
words through a declared depth. It is not a substitute for universal proofs.
"""

from __future__ import annotations

from itertools import product
import sympy as sp

I = sp.eye(2)
Q = sp.Matrix([[1, 1], [1, 0]])
C1 = sp.Matrix([[0, -1], [1, 0]])
Omega = -C1
G = Q.T * Omega - Omega * Q
sqrt5 = sp.sqrt(5)
Phi = (1 + sqrt5) / 2
Psi = (1 - sqrt5) / 2
H = G / sqrt5
S = (2 * Q - I) / sqrt5
v_phi = sp.Matrix([Phi, 1])
v_psi = sp.Matrix([Psi, 1])
P_phi = (Q - Psi * I) / sqrt5
P_psi = (Phi * I - Q) / sqrt5
basis = [I, Q, Omega, G]

checks = 0


def check(condition: bool, label: str) -> None:
    global checks
    if not bool(condition):
        raise AssertionError(label)
    checks += 1


def meq(a: sp.Matrix, b: sp.Matrix) -> bool:
    return all(sp.simplify(x) == 0 for x in (a - b))


# Primitive and quadratic identities
check(meq(C1**2, -I), "C1^2=-I")
check(meq(C1**4, I), "C1^4=I")
check(meq(Omega**2, -I), "Omega^2=-I")
check(meq(Q**2, Q + I), "Q^2=Q+I")
check(meq(Q.inv(), Q - I), "Q^-1=Q-I")
check(meq(G, sp.Matrix([[-2, 1], [1, 2]])), "G value")
check(meq(G, C1 * Q - Q * C1), "G=[C1,Q]")

# Fibonacci powers
F = [sp.Integer(0), sp.Integer(1)]
for _ in range(2, 25):
    F.append(F[-1] + F[-2])
for n in range(1, 24):
    expected = F[n] * Q + F[n - 1] * I
    check(meq(Q**n, expected), f"Fibonacci power n={n}")

# Basis determinant and normal-form extraction
B = sp.Matrix.hstack(
    sp.Matrix([1, 0, 0, 1]),
    sp.Matrix([1, 1, 1, 0]),
    sp.Matrix([0, 1, -1, 0]),
    sp.Matrix([-2, 1, 1, 2]),
)
check(B.det() == -10, "basis determinant")

a, b, c, d = sp.symbols("a b c d")
x0 = (2*a + 3*d - b - c) / 5
x1 = (a - d + 2*b + 2*c) / 5
x2 = (b - c) / 2
x3 = (-2*a + 2*d + b + c) / 10
A_symbolic = sp.Matrix([[a, b], [c, d]])
A_nf = x0*I + x1*Q + x2*Omega + x3*G
check(meq(A_symbolic, A_nf), "symbolic normal-form extraction")

# Multiplication table
relations = [
    (Q*Q, I + Q, "Q^2"),
    (Omega*Omega, -I, "Omega^2"),
    (G*G, 5*I, "G^2"),
    (Q*Omega, (Omega + G)/2, "QOmega"),
    (Omega*Q, (Omega - G)/2, "OmegaQ"),
    (Q*G, (G + 5*Omega)/2, "QG"),
    (G*Q, (G - 5*Omega)/2, "GQ"),
    (Omega*G, 2*Q - I, "OmegaG"),
    (G*Omega, I - 2*Q, "GOmega"),
]
for lhs, rhs, label in relations:
    check(meq(lhs, rhs), label)

# Symbolic coordinate multiplication
x = sp.symbols("x0:4")
y = sp.symbols("y0:4")
A = sum((x[i] * basis[i] for i in range(4)), sp.zeros(2))
Bmat = sum((y[i] * basis[i] for i in range(4)), sp.zeros(2))
z0 = x[0]*y[0] + x[1]*y[1] - x[2]*y[2] - x[2]*y[3] + x[3]*y[2] + 5*x[3]*y[3]
z1 = x[0]*y[1] + x[1]*y[0] + x[1]*y[1] + 2*x[2]*y[3] - 2*x[3]*y[2]
z2 = x[0]*y[2] + x[2]*y[0] + (x[1]*y[2] + x[2]*y[1])/2 + 5*(x[1]*y[3] - x[3]*y[1])/2
z3 = x[0]*y[3] + x[3]*y[0] + (x[1]*y[2] + x[1]*y[3] - x[2]*y[1] + x[3]*y[1])/2
AB_nf = z0*I + z1*Q + z2*Omega + z3*G
check(meq(A*Bmat, AB_nf), "symbolic coordinate multiplication")

# Spectral/null identities
check(meq(Q*v_phi, Phi*v_phi), "Phi eigenvector")
check(meq(Q*v_psi, Psi*v_psi), "Psi eigenvector")
check(sp.simplify(Phi + Psi - 1) == 0, "Phi+Psi=1")
check(sp.simplify(Phi*Psi + 1) == 0, "PhiPsi=-1")
check(sp.simplify(Phi - Psi - sqrt5) == 0, "Phi-Psi=sqrt5")
check(sp.simplify((v_phi.T * Omega * v_psi)[0] - sqrt5) == 0, "oriented separation")
xx, yy = sp.symbols("x y")
v = sp.Matrix([xx, yy])
delta = sp.expand((v.T * G * v)[0])
factored = 2*(yy - xx/Phi)*(yy + Phi*xx)
check(sp.simplify(delta - factored) == 0, "Delta_G factorization")

# Derived projective reflection M0
M0 = Q * Omega
check(meq(M0, sp.Matrix([[-1, 1], [0, 1]])), "M0 value")
check(meq(M0**2, I), "M0 involution")
check(meq(M0*C1, Q), "Q=M0 C1")

# Reflections, Clifford, and group relations
check(meq(G**2, 5*I), "G^2=5I")
check(meq(H**2, I), "H^2=I")
check(meq(H.T*H, I), "H orthogonal")
check(sp.simplify(H.det() + 1) == 0, "det H=-1")
check(meq(S**2, I), "S^2=I")
check(meq(S.T*S, I), "S orthogonal")
check(sp.simplify(S.det() + 1) == 0, "det S=-1")
check(meq(H*Omega, -Omega*H), "H anticommutes Omega")
check(meq(S*H, Omega), "SH=Omega")
check(meq(H*S, -Omega), "HS=-Omega")
check(meq(H*Omega*H, Omega.inv()), "H Omega H=Omega^-1")
check(meq(H*Q*H, -Q.inv()), "H Q H=-Q^-1")
check(meq(H*v_phi, Phi*v_psi), "H exchanges Phi ray")
check(meq(H*v_psi, (1/Phi)*v_phi), "H exchanges Psi ray")
check(meq(Q, (I + Omega*G)/2), "recover Q")

# Projectors and normalized recurrence
check(meq(P_phi + P_psi, I), "projectors sum")
check(meq(P_phi*P_psi, sp.zeros(2)), "projectors orthogonal")
check(meq(P_phi**2, P_phi), "P_phi idempotent")
check(meq(P_psi**2, P_psi), "P_psi idempotent")
check(meq(Q*P_phi, Phi*P_phi), "P_phi spectral")
check(meq(Q*P_psi, Psi*P_psi), "P_psi spectral")
for n in range(0, 16):
    lhs = Phi**(-n) * Q**n
    rhs = P_phi + (-Phi**(-2))**n * P_psi
    check(meq(lhs, rhs), f"normalized recurrence n={n}")

# Exact reduction of all words through depth 9
alphabet = {"Q": Q, "C": C1}
for length in range(0, 10):
    for letters in product(alphabet, repeat=length):
        W = I
        for letter in letters:
            W = W * alphabet[letter]
        aa, bb, cc, dd = W[0, 0], W[0, 1], W[1, 0], W[1, 1]
        coords = [
            (2*aa + 3*dd - bb - cc)/5,
            (aa - dd + 2*bb + 2*cc)/5,
            (bb - cc)/2,
            (-2*aa + 2*dd + bb + cc)/10,
        ]
        reconstructed = sum((coords[i]*basis[i] for i in range(4)), sp.zeros(2))
        check(meq(W, reconstructed), f"word reduction {''.join(letters) or 'empty'}")

print(f"PASS: {checks} exact and finite-regression checks")
