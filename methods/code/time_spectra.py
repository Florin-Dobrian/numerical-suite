"""
TIME_INTEGRATION.md sec.1.1 -- where the modes come from.

For each model operator the eigenvalues are available in closed form. This script
checks each formula against the assembled matrix, tabulates how the spectrum
scales with the mesh, and draws the picture that decides stability:
the scaled spectrum dt*lambda must lie inside the scheme's stability region.

Notation, matching the document:
    M      number of modes = number of unknowns: n-1 for Dirichlet on n cells,
           n for periodic. The +-1 is why |lambda|_max falls just short of its
           limiting value 4*alpha/dx^2 on the Dirichlet problems.
    m      integer mode number, the ONLY index: m = 1 .. M
    k_m    wavenumber of mode m, determined by m: m*pi (Dirichlet), 2*pi*m (periodic)
    phi_m  = k_m*dx, the phase shift per cell, 0 < phi <= pi (pi = the sawtooth).
           Named phi, not theta: theta is the theta-method parameter.
    The wave system is 2M x 2M and gives a +-i*omega PAIR per spatial shape.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLU, ORA, GRN, RED, PUR, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#805ad5", "#718096"
np.set_printoptions(precision=6, suppress=True)


# =====================================================================
# 1. closed-form spectra, each checked against the assembled matrix
# =====================================================================
def spec_diffusion(n, alpha=1.0):
    """Dirichlet, interior nodes. lambda_k = -(4a/dx^2) sin^2(k pi dx / 2)."""
    dx = 1.0 / n
    m = n - 1
    k = np.arange(1, m + 1)
    analytic = -(4 * alpha / dx**2) * np.sin(k * np.pi * dx / 2) ** 2
    L = (alpha / dx**2) * (np.diag(-2 * np.ones(m)) + np.diag(np.ones(m - 1), 1)
                           + np.diag(np.ones(m - 1), -1))
    return np.sort(analytic), np.sort(np.linalg.eigvals(L).real)


def spec_advect_central(n, c=1.0):
    """Periodic, central. lambda_k = -i (c/dx) sin(k dx)."""
    dx = 1.0 / n
    th = 2 * np.pi * np.arange(n) / n
    analytic = -1j * (c / dx) * np.sin(th)
    L = np.zeros((n, n))
    for j in range(n):
        L[j, (j + 1) % n] = -c / (2 * dx)
        L[j, (j - 1) % n] = c / (2 * dx)
    return analytic, np.linalg.eigvals(L)


def spec_advect_upwind(n, c=1.0):
    """Periodic, first-order upwind. lambda_k = -(c/dx)(1 - e^{-i k dx})."""
    dx = 1.0 / n
    th = 2 * np.pi * np.arange(n) / n
    analytic = -(c / dx) * (1 - np.exp(-1j * th))
    L = np.zeros((n, n))
    for j in range(n):
        L[j, j] = -c / dx
        L[j, (j - 1) % n] = c / dx
    return analytic, np.linalg.eigvals(L)


def spec_wave(n, c=1.0):
    """u_tt = c^2 u_xx as a first-order system. lambda = +- i w_k,
       w_k = (2c/dx) sin(k pi dx / 2)."""
    dx = 1.0 / n
    m = n - 1
    k = np.arange(1, m + 1)
    w = (2 * c / dx) * np.sin(k * np.pi * dx / 2)
    analytic = np.concatenate([1j * w, -1j * w])
    K = (c**2 / dx**2) * (np.diag(-2 * np.ones(m)) + np.diag(np.ones(m - 1), 1)
                          + np.diag(np.ones(m - 1), -1))
    A = np.block([[np.zeros((m, m)), np.eye(m)], [K, np.zeros((m, m))]])
    return analytic, np.linalg.eigvals(A)


print("=" * 78)
print("1.1a  Closed-form spectra, checked against the assembled operator (n = 32)")
print("=" * 78)
n = 32
for name, fn, formula in (
        ("diffusion (Dirichlet)", spec_diffusion,
         "lambda_k = -(4 alpha/dx^2) sin^2(k pi dx/2)      real, negative"),
        ("advection, central", spec_advect_central,
         "lambda_k = -i (c/dx) sin(k dx)                   PURELY IMAGINARY"),
        ("advection, upwind", spec_advect_upwind,
         "lambda_k = -(c/dx)(1 - e^{-i k dx})              circle, Re < 0"),
        ("wave equation", spec_wave,
         "lambda_k = +- i (2c/dx) sin(k pi dx/2)           PURELY IMAGINARY, in pairs")):
    a, num = fn(n)
    # conjugate pairs share a real part, and roundoff in the numerical real part
    # would reorder them under a plain lexicographic sort -- round before sorting
    key = lambda v: np.array(sorted(v, key=lambda w: (round(w.real, 8), round(w.imag, 8))))
    err = np.max(np.abs(key(a) - key(num)))
    print(f"\n{name}")
    print(f"   {formula}")
    print(f"   max |analytic - numerical| = {err:.2e}")
    print(f"   |lambda|_max = {np.max(np.abs(a)):.2f}   "
          f"max Re = {np.max(a.real):.3f}   max |Im| = {np.max(abs(a.imag)):.2f}")
print("\nreaction  R(u) = -K u:  lambda = -K, a single point, INDEPENDENT of dx")

# =====================================================================
# 2. how the spectrum scales with the mesh
# =====================================================================
print()
print("=" * 78)
print("1.1b  |lambda|_max against mesh size — the exponent is what matters")
print("=" * 78)
print(f"{'n':>6} {'dx':>9} {'diffusion':>13} {'ratio':>7} {'advection':>13} {'ratio':>7}"
      f" {'wave':>11} {'ratio':>7} {'reaction':>10}")
pd = pa = pw = None
for n in (16, 32, 64, 128, 256):
    dx = 1.0 / n
    ld = np.max(np.abs(spec_diffusion(n)[0]))
    la = np.max(np.abs(spec_advect_upwind(n)[0]))
    lw = np.max(np.abs(spec_wave(n)[0]))
    rd = "" if pd is None else f"{ld/pd:.2f}"
    ra = "" if pa is None else f"{la/pa:.2f}"
    rw = "" if pw is None else f"{lw/pw:.2f}"
    print(f"{n:>6} {dx:>9.5f} {ld:>13.1f} {rd:>7} {la:>13.1f} {ra:>7} {lw:>11.1f} {rw:>7}"
          f" {10.0:>10.1f}")
    pd, pa, pw = ld, la, lw
print("\ndiffusion quadruples when the mesh doubles -> |lambda| ~ 1/dx^2")
print("advection and the wave equation merely double     -> |lambda| ~ 1/dx")
print("reaction does not move at all                     -> |lambda| ~ 1")
print("\nTHIS is why diffusion becomes stiff under refinement and advection does not.")

# =====================================================================
# 2b. why the eigenvalues sit where they do: operator symmetry
# =====================================================================
print()
print("=" * 78)
print("1.1d  Symmetry of L decides the DIRECTION of the spectrum")
print("=" * 78)


def periodic_ops(n=8, dx=None, c=1.0, alpha=1.0):
    dx = dx if dx else 1.0 / n
    D = np.zeros((n, n)); C = np.zeros((n, n)); U = np.zeros((n, n))
    for j in range(n):
        D[j, j] = -2 * alpha / dx**2
        D[j, (j + 1) % n] = alpha / dx**2
        D[j, (j - 1) % n] = alpha / dx**2
        C[j, (j + 1) % n] = -c / (2 * dx)
        C[j, (j - 1) % n] = c / (2 * dx)
        U[j, j] = -c / dx
        U[j, (j - 1) % n] = c / dx
    return D, C, U


D, C, U = periodic_ops()
print(f"{'operator':<24} {'L = L^T':>9} {'L = -L^T':>10} {'normal':>8} "
      f"{'eigenvalues are':>20}")
for nm, L in (("diffusion (periodic)", D), ("advection, central", C),
              ("advection, upwind", U)):
    sym = np.allclose(L, L.T)
    skew = np.allclose(L, -L.T)
    normal = np.allclose(L @ L.T, L.T @ L)
    ev = np.linalg.eigvals(L)
    kind = ("REAL" if np.allclose(ev.imag, 0) else
            "PURELY IMAGINARY" if np.allclose(ev.real, 0) else "complex")
    print(f"{nm:<24} {str(sym):>9} {str(skew):>10} {str(normal):>8} {kind:>20}")
print("\nsymmetric  -> real spectrum        (self-adjoint operator: d2/dx2)")
print("antisymmetric -> imaginary spectrum  (skew-adjoint operator: d/dx, periodic)")
print("neither    -> complex spectrum       (upwind: a central part plus a symmetric part)")
print("\nAll three are NORMAL here, which is why the mode picture is exact for them.")
print("Inflow/outflow boundaries and variable coefficients break normality.")

# =====================================================================
# 3. the picture: dt * spectrum must lie inside the stability region
# =====================================================================
fig, axes = plt.subplots(1, 4, figsize=(14, 4.1))
x = np.linspace(-3.2, 1.2, 500)
y = np.linspace(-2.6, 2.6, 500)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y
stable_fe = np.abs(1 + Z) <= 1

n = 32
cases = [
    ("diffusion\n$r = \\alpha\\Delta t/\\Delta x^2 = 0.5$", spec_diffusion(n)[0],
     0.5 / (1.0 * n**2), BLU, "STABLE\nfills the disc exactly,\nreaching $-2$"),
    ("advection, central\n$\\nu = c\\Delta t/\\Delta x = 0.3$", spec_advect_central(n)[0],
     0.3 / n, RED, "NEVER STABLE\nthe disc has zero width\non this axis"),
    ("advection, upwind\n$\\nu = 0.8$", spec_advect_upwind(n)[0], 0.8 / n, GRN,
     "STABLE\na circle inside\nthe disc"),
    ("wave equation\n$\\nu = 0.5$", spec_wave(n)[0], 0.5 / n, PUR,
     "NEVER STABLE\nsame axis, same\nreason"),
]
for i, (ax, (title, lam, dt, col, note)) in enumerate(zip(axes, cases)):
    ax.contourf(X, Y, stable_fe.astype(float), levels=[0.5, 1.5],
                colors=["#cbd5e0"], alpha=0.85)
    ax.contour(X, Y, np.abs(1 + Z), levels=[1.0], colors=["#4a5568"], linewidths=1.2)
    z = lam * dt
    ax.plot(z.real, z.imag, "o", ms=3.5, color=col, zorder=5)
    ax.axhline(0, color=GRY, lw=0.5); ax.axvline(0, color=GRY, lw=0.5)
    ax.set_xlim(-3.3, 1.3); ax.set_ylim(-2.9, 2.9); ax.set_aspect("equal")
    ax.set_title(title, fontsize=8.5)
    ax.text(-1.0, -2.05, note, fontsize=7, color=col, ha="center", va="center",
            linespacing=1.35)
    ax.set_xlabel(r"Re$(z) = $Re$(\lambda\Delta t)$", fontsize=7.5)
    if i == 0:
        ax.annotate("grey disc:\nforward Euler's\nstability region\n$|1+z|\\leq 1$",
                    xy=(-1.0, 0.85), xytext=(-3.2, 2.75), fontsize=7, color="#2d3748",
                    ha="left", va="top", linespacing=1.3,
                    arrowprops=dict(arrowstyle="->", color="#2d3748", lw=0.8))
        ax.annotate("coloured dots:\n$z_m=\\lambda_m\\Delta t$,\none per mode",
                    xy=(-1.75, 0.0), xytext=(-3.2, -0.75), fontsize=7, color=BLU,
                    ha="left", va="top", linespacing=1.3,
                    arrowprops=dict(arrowstyle="->", color=BLU, lw=0.8))
axes[0].set_ylabel(r"Im$(z) = $Im$(\lambda\Delta t)$", fontsize=7.5)
fig.suptitle("Forward Euler is stable if and only if EVERY scaled eigenvalue "
             "$z_m=\\lambda_m\\Delta t$ lies inside its stability region",
             y=1.04, fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/time_fig9_spectra.png", bbox_inches="tight")
plt.close(fig)

print()
print("=" * 78)
print("1.1c  Does dt * spectrum fit inside the forward Euler region?  (n = 32)")
print("=" * 78)
print(f"{'operator':<24} {'parameter':>14} {'max |1 + z|':>13}   verdict")
for label, lam, dt, par in (
        ("diffusion", spec_diffusion(n)[0], 0.5 / n**2, "r = 0.50"),
        ("diffusion", spec_diffusion(n)[0], 0.51 / n**2, "r = 0.51"),
        ("advection, central", spec_advect_central(n)[0], 0.3 / n, "nu = 0.30"),
        ("advection, central", spec_advect_central(n)[0], 0.01 / n, "nu = 0.01"),
        ("advection, upwind", spec_advect_upwind(n)[0], 0.8 / n, "nu = 0.80"),
        ("advection, upwind", spec_advect_upwind(n)[0], 1.05 / n, "nu = 1.05"),
        ("wave equation", spec_wave(n)[0], 0.5 / n, "nu = 0.50")):
    g = np.max(np.abs(1 + lam * dt))
    print(f"{label:<24} {par:>14} {g:>13.6f}   "
          f"{'UNSTABLE' if g > 1 + 1e-12 else 'stable'}")
print("\nCentral advection fails at nu = 0.01 as surely as at nu = 0.30: its modes sit ON")
print("the imaginary axis, and the forward Euler region touches that axis only at zero.")
print(f"\nfigure -> {os.path.normpath(FIG)}/time_fig9_spectra.png")
