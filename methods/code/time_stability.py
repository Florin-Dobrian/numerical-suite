"""
TIME_INTEGRATION.md sec.1 -- stiffness without space.

The Dahlquist test equation y' = lambda y isolates every stability question
from every spatial question. Run:  uv run time_stability.py
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.set_printoptions(precision=6, suppress=True, linewidth=140)
BLU, ORA, GRN, RED, PUR, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#805ad5", "#718096"

# --- amplification factors R(z), z = lambda*dt -------------------------------
R = {
    "Forward Euler":  lambda z: 1 + z,
    "Backward Euler": lambda z: 1 / (1 - z),
    "Crank-Nicolson": lambda z: (1 + z / 2) / (1 - z / 2),
    "RK4":            lambda z: 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24,
}
ORDER = {"Forward Euler": 1, "Backward Euler": 1, "Crank-Nicolson": 2, "RK4": 4}

print("=" * 78)
print("1.1  Amplification factor and the two kinds of stability")
print("=" * 78)
print(f"{'scheme':<16} {'order':>5} {'R(-0.5)':>10} {'R(-2)':>10} {'R(-1000)':>12} "
      f"{'A-stable':>9} {'L-stable':>9}")
for name, f in R.items():
    lhp = np.array([f(z) for z in
                    (-0.5 + 0j, -2 + 0j, -1e3 + 0j, -1e6 + 0j,
                     -1 + 10j, -0.1 + 100j, -1e3 + 1e3j)])
    A = bool(np.all(np.abs(lhp) <= 1 + 1e-12))
    Linf = abs(f(-1e12 + 0j))
    print(f"{name:<16} {ORDER[name]:>5} {f(-0.5+0j).real:>10.4f} {f(-2+0j).real:>10.4f} "
          f"{f(-1e3+0j).real:>12.6f} {str(A):>9} {str(Linf < 1e-8):>9}")

print("\nThe L-stability failure of Crank-Nicolson, concretely.")
print("A mode with lambda*dt = -1000 -- fast, and physically dead within one step:")
for name in ("Backward Euler", "Crank-Nicolson"):
    r = R[name](-1000 + 0j).real
    print(f"  {name:<16} R = {r:+.6f}   after 50 steps: {r**50:+.6e}")
print("  Backward Euler annihilates it. Crank-Nicolson flips its sign and keeps")
print("  almost all of it -- a stable, oscillating, entirely spurious mode.")

# --- 1.2 stability regions ---------------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(13, 3.4))
x = np.linspace(-6, 4, 700)
y = np.linspace(-5, 5, 700)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y
for a, (name, f) in zip(axes, R.items()):
    with np.errstate(all="ignore"):
        M = np.abs(f(Z))
    a.contourf(X, Y, (M <= 1).astype(float), levels=[0.5, 1.5], colors=[BLU], alpha=0.35)
    a.contour(X, Y, M, levels=[1.0], colors=[BLU], linewidths=1.4)
    a.axhline(0, color=GRY, lw=0.6); a.axvline(0, color=GRY, lw=0.6)
    a.axvspan(-6, 0, color=GRN, alpha=0.06)
    a.set_xlim(-6, 4); a.set_ylim(-5, 5); a.set_aspect("equal")
    a.set_title(f"{name}\norder {ORDER[name]}", fontsize=9)
    a.set_xlabel(r"Re$(\lambda\Delta t)$", fontsize=8)
axes[0].set_ylabel(r"Im$(\lambda\Delta t)$", fontsize=8)
fig.suptitle("Stability regions  $|R(z)| \\leq 1$  (blue).  Green: the left half plane, "
             "where the true solution decays.", y=1.04, fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/time_fig1_stability_regions.png", bbox_inches="tight")
plt.close(fig)

# --- 1.3 order verification on a non-stiff problem ---------------------------
print()
print("=" * 78)
print("1.3  Observed order of accuracy,  y' = -y,  y(0)=1,  to t=1")
print("=" * 78)
lam, T = -1.0, 1.0
exact = np.exp(lam * T)


def integrate(name, dt):
    n = int(round(T / dt)); y = 1.0
    if name == "RK4":
        for _ in range(n):
            k1 = lam * y; k2 = lam * (y + dt * k1 / 2)
            k3 = lam * (y + dt * k2 / 2); k4 = lam * (y + dt * k3)
            y += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    else:
        g = R[name](lam * dt).real
        y = g**n
    return y


dts = [0.2, 0.1, 0.05, 0.025, 0.0125]
print(f"{'dt':>8} " + "".join(f"{k:>26}" for k in R))
prev = {k: None for k in R}
for dt in dts:
    row = f"{dt:>8.4f} "
    for k in R:
        e = abs(integrate(k, dt) - exact)
        o = "" if prev[k] is None else f" (p={np.log(prev[k]/e)/np.log(2):.2f})"
        row += f"{e:>14.3e}{o:>12}"
        prev[k] = e
    print(row)

# --- 1.4 a genuinely stiff system -------------------------------------------
print()
print("=" * 78)
print("1.4  A stiff system: eigenvalues -1 and -1000, integrated to t = 5")
print("=" * 78)
A = np.array([[-1.0, 0.0], [999.0, -1000.0]])      # eigenvalues -1, -1000
y0 = np.array([1.0, 0.0])
Tend = 5.0
evals = np.linalg.eigvals(A)
print(f"eigenvalues: {evals}      stiffness ratio S = {abs(evals).max()/abs(evals).min():.0f}")
print(f"slow timescale 1/|lambda_slow| = {1/abs(evals).min():.3f} s -- this is what we care about")
print(f"fast timescale 1/|lambda_fast| = {1/abs(evals).max():.5f} s -- dead after ~0.005 s\n")

dt_fe_max = 2.0 / abs(evals).max()
print(f"Forward Euler stability limit dt <= 2/|lambda_max| = {dt_fe_max:.6f} s")
print(f"  -> steps required to reach t = {Tend}:  {int(np.ceil(Tend/dt_fe_max)):,}")
print(f"Backward Euler at dt = 0.05 (accuracy-limited, not stability-limited)")
print(f"  -> steps required:  {int(Tend/0.05):,}    ratio {int(np.ceil(Tend/dt_fe_max))/(Tend/0.05):.0f}x fewer\n")

I = np.eye(2)
for label, dt in (("dt = 0.0019 (just inside the limit)", 0.0019),
                  ("dt = 0.0021 (just outside)", 0.0021)):
    y = y0.copy()
    for _ in range(int(Tend / dt)):
        y = y + dt * (A @ y)
        if not np.all(np.isfinite(y)) or np.abs(y).max() > 1e10:
            break
    print(f"  Forward Euler, {label:<36} |y| = {np.abs(y).max():.3e} "
          f"{'DIVERGED' if np.abs(y).max() > 1e10 else 'ok'}")
for dt in (0.05, 0.5):
    y = y0.copy(); M = np.linalg.inv(I - dt * A)
    for _ in range(int(Tend / dt)):
        y = M @ y
    print(f"  Backward Euler, dt = {dt:<41} |y| = {np.abs(y).max():.3e} ok")

print("\nThe fast mode is stable, physically irrelevant after 5 ms, and it alone")
print("sets the explicit step for all 5 seconds. That is stiffness.")
print(f"\nfigures -> {os.path.normpath(FIG)}")
