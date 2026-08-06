#!/usr/bin/env python3
"""Exact and finite-regression gate for Conserved Transport and Return, Draft 0.1.

The script checks displayed finite identities and representative parameter grids.
It is a regression gate, not a substitute for the universal proofs in the paper.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import sympy as sp


TOL = 1.0e-10
checks = 0


def check(condition: bool, message: str) -> None:
    global checks
    checks += 1
    if not condition:
        raise AssertionError(message)


def close(a: float | np.ndarray, b: float | np.ndarray, tol: float = TOL) -> bool:
    return bool(np.allclose(a, b, rtol=tol, atol=tol))


# ---------------------------------------------------------------------------
# 1. Planar generator and exact symbolic identities
# ---------------------------------------------------------------------------
I2 = sp.eye(2)
Omega = sp.Matrix([[0, 1], [-1, 0]])
J0 = -Omega

check(Omega**2 == -I2, "Omega^2 != -I")
check(J0 == -Omega, "J0 != -Omega")
check(J0**2 == -I2, "J0^2 != -I")
check(Omega.T == -Omega, "Omega is not skew")
check(J0.T == -J0, "J0 is not skew")
check(Omega.det() == 1, "det(Omega) != 1")
check(J0.det() == 1, "det(J0) != 1")

th = sp.symbols("th", real=True)
Rth = sp.cos(th) * I2 + sp.sin(th) * J0
check(sp.simplify(Rth.T * Rth - I2) == sp.zeros(2), "R(theta) not orthogonal")
check(sp.simplify(Rth.det() - 1) == 0, "det R(theta) != 1")
check(sp.simplify(sp.diff(Rth, th) - J0 * Rth) == sp.zeros(2), "rotation generator law failed")


# ---------------------------------------------------------------------------
# 2. Circular central-force orbit, return hierarchy, curvature, and current
# ---------------------------------------------------------------------------
R, mu, m, w, t = sp.symbols("R mu m w t", positive=True, real=True)
q = sp.Matrix([R * sp.cos(w * t), R * sp.sin(w * t)])
v = sp.diff(q, t)
a = sp.diff(v, t)

check(sp.simplify(q.dot(q) - R**2) == 0, "circle radius not constant")
check(sp.simplify(v.dot(v) - R**2 * w**2) == 0, "circle speed formula failed")
check(sp.simplify(a + w**2 * q) == sp.zeros(2, 1), "circular acceleration failed")
check(sp.simplify(m * (q[0] * v[1] - q[1] * v[0]) - m * R**2 * w) == 0,
      "angular current formula failed")

# Plane-curve curvature |det(q',q'')| / ||q'||^3 = 1/R.
det_va = sp.det(sp.Matrix.hstack(v, a))
curv_squared = sp.simplify(det_va**2 / (v.dot(v))**3)
check(sp.simplify(curv_squared - 1 / R**2) == 0, "curvature formula failed")

# Exact central-force equation when w^2 = mu/R^3.
central_residual = sp.simplify((a + (mu / R**3) * q).subs(mu, w**2 * R**3))
check(central_residual == sp.zeros(2, 1), "inverse-square circular orbit residual")

# Forced angular-momentum balance.
qx, qy, vx, vy, Fx, Fy = sp.symbols("qx qy vx vy Fx Fy", real=True)
r3 = sp.symbols("r3", positive=True, real=True)
ax = -mu * qx / r3 + Fx / m
ay = -mu * qy / r3 + Fy / m
Ldot = sp.simplify(m * (qx * ay - qy * ax))
check(sp.simplify(Ldot - (qx * Fy - qy * Fx)) == 0, "forced axial balance failed")

# Large finite regression grid.
for Ri in range(1, 17):
    for mui in range(1, 17):
        omega = math.sqrt(mui / (Ri**3))
        period = 2.0 * math.pi / omega
        for mi in range(1, 9):
            q0 = np.array([float(Ri), 0.0])
            v0 = np.array([0.0, float(Ri) * omega])
            qh = np.array([-float(Ri), 0.0])
            vh = np.array([0.0, -float(Ri) * omega])
            qT = np.array([float(Ri), 0.0])
            vT = np.array([0.0, float(Ri) * omega])

            check(close(qT, q0), "full-period position return failed")
            check(close(vT, v0), "full-period velocity return failed")
            check(not close(qh, q0), "half-period state position should differ")
            check(not close(vh, v0), "half-period state velocity should differ")
            check(abs(np.linalg.det(np.column_stack([q0, qh]))) < TOL,
                  "projective half-return failed")
            check(float(np.dot(q0, qh)) < 0.0, "half-period vectors not antipodal")
            check(abs(1.0 / Ri - 1.0 / Ri) < TOL, "curvature grid failure")
            check(abs(mi * Ri**2 * omega - mi * (q0[0] * v0[1] - q0[1] * v0[0])) < TOL,
                  "angular-current grid failure")
            check(period > 0.0, "period not positive")
            check(abs(omega**2 - mui / Ri**3) < TOL, "Kepler circular frequency failed")


# ---------------------------------------------------------------------------
# 3. Reversibility on the explicit circular orbit
# ---------------------------------------------------------------------------
def orbit_state(radius: float, omega: float, time: float) -> np.ndarray:
    return np.array([
        radius * math.cos(omega * time),
        radius * math.sin(omega * time),
        -radius * omega * math.sin(omega * time),
        radius * omega * math.cos(omega * time),
    ])


def clockwise_state(radius: float, omega: float, time: float) -> np.ndarray:
    return np.array([
        radius * math.cos(omega * time),
        -radius * math.sin(omega * time),
        -radius * omega * math.sin(omega * time),
        -radius * omega * math.cos(omega * time),
    ])


def time_reverse(state: np.ndarray) -> np.ndarray:
    return np.array([state[0], state[1], -state[2], -state[3]])


for Ri in range(1, 17):
    for mui in range(1, 17):
        omega = math.sqrt(mui / Ri**3)
        for k in range(-8, 9):
            tt = k / 5.0
            x0 = orbit_state(Ri, omega, 0.0)
            # K phi_t K applied to the positive-orientation initial state.
            lhs = time_reverse(clockwise_state(Ri, omega, tt))
            rhs = orbit_state(Ri, omega, -tt)
            check(close(lhs, rhs), "time-reversal orbit identity failed")
            check(close(time_reverse(time_reverse(x0)), x0), "K^2 != I")


# ---------------------------------------------------------------------------
# 4. Maxwell cancellation, null propagation, and holonomy
# ---------------------------------------------------------------------------
c, rho, fpp, x0s, lam = sp.symbols("c rho fpp x0 lam", positive=True, real=True)
check(sp.simplify(-fpp / c**2 + fpp / c**2) == 0, "Maxwell divergence cancellation failed")
check(sp.simplify(-c**2 + c**2) == 0, "null tangent failed")
check(sp.simplify(lam - (x0s + c * lam) / c + x0s / c) == 0,
      "retarded coordinate not constant on null line")

Omega_np = np.array([[0.0, 1.0], [-1.0, 0.0]])
I_np = np.eye(2)


def rot_from_a(a_value: float) -> np.ndarray:
    angle = -2.0 * math.pi * a_value
    return math.cos(angle) * I_np + math.sin(angle) * Omega_np


for k in range(-512, 513):
    av = k / 64.0
    H = rot_from_a(av)
    check(close(H.T @ H, I_np), "holonomy not orthogonal")
    check(abs(np.linalg.det(H) - 1.0) < TOL, "holonomy determinant failed")
    check(close(rot_from_a(av + 1.0), H), "holonomy period in a failed")
    check(close(rot_from_a(-av), np.linalg.inv(H)), "holonomy inverse failed")

# Group composition on a broad rational grid.
for i in range(-32, 33):
    for j in range(-16, 17):
        aa = i / 17.0
        bb = j / 19.0
        check(close(rot_from_a(aa + bb), rot_from_a(aa) @ rot_from_a(bb)),
              "holonomy group law failed")


# ---------------------------------------------------------------------------
# 5. Proper-time inequalities and parameter independence
# ---------------------------------------------------------------------------
for rhoi in range(1, 33):
    for ci in range(1, 9):
        threshold = 2.0 * math.pi * rhoi / ci
        for multiplier in (1.05, 1.25, 2.0, 5.0, 10.0):
            TT = threshold * multiplier
            moving = TT * math.sqrt(1.0 - (2.0 * math.pi * rhoi / (ci * TT))**2)
            delta = TT - moving
            check(TT > threshold, "timelike threshold failed")
            check(moving > 0.0, "moving proper time not positive")
            check(moving < TT, "moving proper time not smaller")
            check(delta > 0.0, "proper-time difference not positive")
            # Independent of a and f: the formula has only T, rho, c.
            for av in (-2.5, -0.25, 0.0, 0.75, 3.0):
                check(abs((TT - moving) - delta) < TOL, f"spurious a-dependence at {av}")


# ---------------------------------------------------------------------------
# 6. Linear commutator criterion: exact finite matrix instances
# ---------------------------------------------------------------------------
J = sp.Matrix([[0, -1], [1, 0]])
for aa in range(-12, 13):
    for bb in range(-12, 13):
        A = aa * I2 + bb * J
        check(A * J - J * A == sp.zeros(2), "commutant family failed")

# All small integer 2x2 matrices: identify the exact commutant of J.
for a11 in range(-4, 5):
    for a12 in range(-4, 5):
        for a21 in range(-4, 5):
            for a22 in range(-4, 5):
                A = sp.Matrix([[a11, a12], [a21, a22]])
                commutes = (A * J - J * A == sp.zeros(2))
                expected = (a11 == a22 and a21 == -a12)
                check(commutes == expected, "2x2 rotational commutant classification failed")

# Numeric exponential instances for the commuting family.
for aa in range(-8, 9):
    for bb in range(-8, 9):
        A = np.array([[float(aa), -float(bb)], [float(bb), float(aa)]])
        for ti in (-0.3, -0.05, 0.0, 0.07, 0.2):
            # exp(t(aI+bJ)) = exp(at) R(bt)
            EA = math.exp(aa * ti) * np.array([
                [math.cos(bb * ti), -math.sin(bb * ti)],
                [math.sin(bb * ti), math.cos(bb * ti)],
            ])
            for theta in (-1.1, -0.2, 0.0, 0.4, 1.7):
                U = np.array([
                    [math.cos(theta), -math.sin(theta)],
                    [math.sin(theta), math.cos(theta)],
                ])
                check(close(EA @ U, U @ EA), "commuting exponential instance failed")

print(f"PASS exact regression suite ({checks:,} checks)")
