"""
TIME_INTEGRATION.md sec.1.1 -- what a "mode" is, and what makes one fast.

The semi-discrete operator L has eigenvectors. Any solution is a sum over them,
each decaying at its own rate. This script shows what they look like, how fast
each one dies, and how refining the mesh manufactures stiffness.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLU, ORA, GRN, RED, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096"
ALPHA = 1.0


def laplacian_modes(nx, alpha=ALPHA):
    """Eigenpairs of the Dirichlet discrete Laplacian, slowest first."""
    dx = 1.0 / nx
    m = nx - 1
    L = (alpha / dx**2) * (np.diag(-2 * np.ones(m))
                           + np.diag(np.ones(m - 1), 1)
                           + np.diag(np.ones(m - 1), -1))
    w, V = np.linalg.eigh(L)
    o = np.argsort(-w)
    return w[o], V[:, o], dx


nx = 24
w, V, dx = laplacian_modes(nx)
m = nx - 1

print("=" * 78)
print(f"1.1  Eigenmodes of the discrete Laplacian   (nx = {nx}, alpha = {ALPHA},"
      f" dx = {dx:.4f})")
print("=" * 78)
print("A solution is a sum over these:  u(t) = sum_k c_k exp(lambda_k t) v_k,")
print("and each term evolves entirely independently of the others.\n")
print(f"{'mode k':>7} {'lambda_k':>12} {'decay time 1/|lambda|':>22}   shape")
for k in (1, 2, 5, 12, m - 1, m):
    tau = 1 / abs(w[k - 1])
    shape = ("smoothest — one half-sine" if k == 1 else
             "SAWTOOTH — the wiggliest the mesh can hold" if k == m else
             f"{k} half-waves")
    print(f"{k:>7} {w[k-1]:>12.2f} {tau:>22.3e}   {shape}")

print(f"\nSmooth modes are slow, wiggly modes are fast. The extreme is the sawtooth,")
print(f"at lambda ~ -4 alpha / dx^2 = {-4*ALPHA/dx**2:.1f} (measured {w[-1]:.1f}).")
print(f"\nstiffness ratio |lambda_max| / |lambda_min| = {abs(w[-1])/abs(w[0]):.0f}")

print()
print("=" * 78)
print("1.1b  Refining the mesh manufactures stiffness")
print("=" * 78)
print(f"{'nx':>6} {'|lambda_min| (slow)':>21} {'|lambda_max| (fast)':>21} "
      f"{'ratio':>10} {'explicit dt <= 2/|l_max|':>26}")
for n in (12, 24, 48, 96, 192):
    ww, _, dxn = laplacian_modes(n)
    print(f"{n:>6} {abs(ww[0]):>21.2f} {abs(ww[-1]):>21.1f} "
          f"{abs(ww[-1])/abs(ww[0]):>10.0f} {2/abs(ww[-1]):>26.3e}")
print("\nThe slow mode is physics and barely moves: it is the fundamental of the domain.")
print("The fast mode is a property of the GRID and scales as 1/dx^2. Doubling nx")
print("quadruples the stiffness ratio and quarters the explicit step, while the")
print("answer you are after has not changed at all.")

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
x = np.linspace(0, 1, nx + 1)
for k, c, lbl in ((1, BLU, "k=1  slowest"), (2, GRN, "k=2"),
                  (5, ORA, "k=5"), (m, RED, "k=23  FASTEST")):
    v = np.zeros(nx + 1)
    v[1:-1] = V[:, k - 1] / np.max(np.abs(V[:, k - 1]))
    ax[0].plot(x, v, "-o", ms=2.5, color=c, lw=1.3,
               label=f"{lbl}   $\\lambda$={w[k-1]:.0f}")
ax[0].axhline(0, color=GRY, lw=0.7)
ax[0].set_ylim(-1.1, 1.75)
ax[0].set_xlabel("$x$"); ax[0].set_ylabel("eigenvector")
ax[0].set_title("What the modes look like", fontsize=10)
ax[0].legend(fontsize=7.5, frameon=False, ncol=2, loc="upper center")

t = np.linspace(0, 0.02, 400)
for k, c, lbl in ((1, BLU, "k=1"), (2, GRN, "k=2"), (5, ORA, "k=5"), (m, RED, "k=23")):
    ax[1].semilogy(t, np.exp(w[k - 1] * t), color=c, lw=1.5,
                   label=f"{lbl},  $\\tau$={1/abs(w[k-1]):.1e}")
ax[1].set_ylim(1e-6, 2)
ax[1].set_xlabel("$t$"); ax[1].set_ylabel(r"amplitude $e^{\lambda_k t}$")
ax[1].set_title("How fast each one dies", fontsize=10)
ax[1].legend(fontsize=7.5, frameon=False); ax[1].grid(alpha=0.25)

fig.suptitle("Eigenmodes of the discrete Laplacian — 'fast' means 'wiggly'",
             y=1.03, fontsize=10.5)
fig.tight_layout(); fig.savefig(f"{FIG}/time_fig8_modes.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigure -> {os.path.normpath(FIG)}/time_fig8_modes.png")
