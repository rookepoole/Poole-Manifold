#!/usr/bin/env python3
"""Poole Manifold PM0: no-cost, zero-dependency reference verifier.

Python standard library only. No Mathematica, Wolfram Cloud, pip, GPU,
accounts or network access required.

Checks the synchronous B5-7/S5-9 26-neighbor update against the
independently derived q_D/r_D first-step sheet-defect predictor.
"""
from __future__ import annotations
from random import Random

OFFSETS3 = [(a, b, c) for a in (-1, 0, 1)
                      for b in (-1, 0, 1) for c in (-1, 0, 1)
                      if (a, b, c) != (0, 0, 0)]
OFFSETS2 = [(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)
            if (a, b) != (0, 0)]


def make_sheet(defect):
    L = len(defect)
    assert L >= 5 and len(defect[0]) == L
    x = [[[0 for _ in range(L)] for _ in range(L)] for _ in range(L)]
    z = L // 2
    for j in range(L):
        for k in range(L):
            x[z][j][k] = 1 - defect[j][k]
    return x


def rule_step(x, survival_upper_bonus=2):
    """Literal full 3D update, independent of the planar q/r formulas."""
    L = len(x)
    y = [[[0 for _ in range(L)] for _ in range(L)] for _ in range(L)]
    for i in range(L):
        for j in range(L):
            for k in range(L):
                n = sum(x[(i + di) % L][(j + dj) % L][(k + dk) % L]
                        for di, dj, dk in OFFSETS3)
                y[i][j][k] = int(5 <= n <= 7 + survival_upper_bonus * x[i][j][k])
    return y


def planar_prediction(defect):
    """Analytic first-step prediction from two-dimensional missing-cell counts."""
    L = len(defect)
    y = [[[0 for _ in range(L)] for _ in range(L)] for _ in range(L)]
    z = L // 2
    for j in range(L):
        for k in range(L):
            q = sum(defect[(j + dj) % L][(k + dk) % L] for dj, dk in OFFSETS2)
            r = q + defect[j][k]
            y[z][j][k] = int(q <= 3 if defect[j][k] == 0 else 1 <= q <= 3)
            y[z - 1][j][k] = int(2 <= r <= 4)
            y[z + 1][j][k] = int(2 <= r <= 4)
    return y


def rectangle_mask(L, a, b):
    assert 2 <= a <= L - 10 and 2 <= b <= L - 11
    d = [[0 for _ in range(L)] for _ in range(L)]
    for j in range(5, 5 + a):
        for k in range(6, 6 + b):
            d[j][k] = 1
    return d


def first_step_metrics(x, y):
    L = len(x)
    births = deaths = 0
    channels = {5: 0, 6: 0, 7: 0}
    for i in range(L):
        for j in range(L):
            for k in range(L):
                if x[i][j][k] == 0 and y[i][j][k] == 1:
                    births += 1
                    n = sum(x[(i + di) % L][(j + dj) % L][(k + dk) % L]
                            for di, dj, dk in OFFSETS3)
                    channels[n] = channels.get(n, 0) + 1
                if x[i][j][k] == 1 and y[i][j][k] == 0:
                    deaths += 1
    return births, deaths, tuple(channels[k] for k in (5, 6, 7))


def run():
    assert len(OFFSETS3) == 26
    L = 17
    d = [[0 for _ in range(L)] for _ in range(L)]
    plane = make_sheet(d)
    assert rule_step(plane) == plane
    assert sum(v for layer in rule_step(plane, 0) for row in layer for v in row) == 0
    d[8][8] = 1
    pore = make_sheet(d)
    assert rule_step(pore) == pore
    print('PASS: complete sheet fixed; isolated pore fixed; ablation removes sheet')

    for a, b in [(2, 2), (3, 4), (10, 10), (12, 12)]:
        L = max(25, a + 12, b + 13)
        d = rectangle_mask(L, a, b)
        x = make_sheet(d)
        y = rule_step(x)
        assert y == planar_prediction(d), (a, b, 'voxel mismatch')
        got = first_step_metrics(x, y)
        expected = (4 * a + 4 * b + 12, 0, (12, 4 * a + 4 * b - 16, 16))
        assert got == expected, (a, b, got, expected)
        print(f'PASS: {a}x{b} hole: {got[0]} births, {got[1]} deaths, channels {got[2]}')

    rng = Random(20261018)
    tested = 12
    for case in range(tested):
        density = [0.03, 0.10, 0.20, 0.35, 0.50, 0.70][case % 6]
        d = [[int(rng.random() < density) for _ in range(13)] for _ in range(13)]
        assert rule_step(make_sheet(d)) == planar_prediction(d), ('random mask', case)
    print(f'PASS: {tested} random masks, voxel-for-voxel full 3D vs. 2D predictor')
    print('ALL CHECKS PASSED. No libraries, accounts, or uploads required.')


if __name__ == '__main__':
    run()
