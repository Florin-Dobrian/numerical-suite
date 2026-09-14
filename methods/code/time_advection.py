"""
TIME_INTEGRATION.md sec.3 -- hyperbolic.  u_t + c u_x = 0 on a periodic [0,1], c = 1.

The exact solution translates the initial profile rigidly. Anything a scheme does
other than translate it is error, and it comes in exactly two kinds.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLU, ORA, GRN, RED, PUR, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#805ad5", "#718096"
C = 1.0


def advance(u0, nu, nsteps, scheme, nx):
    u = u0.copy()
    if scheme == "implicit upwind":
        e = np.ones(nx)
        A = sp.diags([1 + nu, -nu * e[:-1], -nu], [0, -1, nx - 1],
                     shape=(nx, nx), format="csc")
        lu = spla.splu(A)
        for _ in range(nsteps):
            u = lu.solve(u)
        return u
    for _ in range(nsteps):
        up, dn = np.roll(u, -1), np.roll(u, 1)          # u_{j+1}, u_{j-1}
        if scheme == "FTCS":
            u = u - 0.5 * nu * (up - dn)
        elif scheme == "upwind":
            u = u - nu * (u - dn)
        elif scheme == "Lax-Friedrichs":
            u = 0.5 * (up + dn) - 0.5 * nu * (up - dn)
        elif scheme == "Lax-Wendroff":
            u = u - 0.5 * nu * (up - dn) + 0.5 * nu**2 * (up - 2 * u + dn)
        else:
            raise ValueError(scheme)
        if not np.all(np.isfinite(u)) or np.abs(u).max() > 1e8:
            return np.full(nx, np.nan)
    return u


nx = 200
x = (np.arange(nx) + 0.5) / nx
gauss = np.exp(-((x - 0.3) / 0.06) ** 2)
square = ((x > 0.2) & (x < 0.4)).astype(float)
SCHEMES = ["FTCS", "upwind", "Lax-Friedrichs", "Lax-Wendroff", "implicit upwind"]

print("=" * 78)
print("3.1  One full revolution of a Gaussian, nx = 200")
print("=" * 78)
print("The exact solution after one period is the initial condition, unchanged.\n")
print(f"{'scheme':<18} {'nu':>6} {'peak':>9} {'L2 error':>11} {'min u':>9}   comment")
res31 = {}
for scheme in SCHEMES:
    for nu in (0.5, 0.9, 1.0, 5.0):
        if scheme != "implicit upwind" and nu > 1.0:
            continue
        nsteps = int(round(1.0 / (nu / nx)))
        u = advance(gauss, nu, nsteps, scheme, nx)
        if not np.all(np.isfinite(u)):
            print(f"{scheme:<18} {nu:>6.2f} {'--':>9} {'--':>11} {'--':>9}   BLEW UP")
            continue
        err = np.sqrt(np.mean((u - gauss) ** 2))
        note = ""
        if scheme == "upwind" and abs(nu - 1) < 1e-12:
            note = "exact: a one-cell shift per step"
        if scheme == "implicit upwind":
            note = "stable, and badly wrong"
        print(f"{scheme:<18} {nu:>6.2f} {u.max():>9.4f} {err:>11.3e} {u.min():>9.4f}   {note}")
        res31[(scheme, nu)] = u

print("\nFTCS is unconditionally unstable: |g|^2 = 1 + nu^2 sin^2(k dx) > 1 for every")
print("nu > 0 and every mode. No time step saves it -- this is not a CFL violation.")

print()
print("=" * 78)
print("3.2  The two kinds of error, on a square pulse (one revolution, nu = 0.5)")
print("=" * 78)
nsteps = int(round(1.0 / (0.5 / nx)))
print(f"{'scheme':<18} {'max u':>9} {'min u':>9} {'undershoot':>12} {'total variation':>17}")
print(f"{'exact':<18} {1.0:>9.4f} {0.0:>9.4f} {0.0:>12.4f} {2.0:>17.4f}")
res32 = {}
for scheme in ("upwind", "Lax-Wendroff", "Lax-Friedrichs"):
    u = advance(square, 0.5, nsteps, scheme, nx)
    res32[scheme] = u
    print(f"{scheme:<18} {u.max():>9.4f} {u.min():>9.4f} {min(0.0, u.min()):>12.4f} "
          f"{np.sum(np.abs(np.diff(np.r_[u, u[0]]))):>17.4f}")
print("\nupwind      -- monotone, no undershoot, and the pulse is smeared: DISSIPATION")
print("Lax-Wendroff-- amplitude kept, oscillations behind the jump:    DISPERSION")
print("These are different failures and neither implies the other.")

print()
print("=" * 78)
print("3.3  Modified equation: upwind is exact advection plus artificial diffusion")
print("=" * 78)
print("Taylor expanding the upwind update gives")
print("    u_t + c u_x = alpha_num u_xx  with  alpha_num = c dx (1 - nu) / 2")
print("so the artificial diffusivity vanishes at nu = 1 and is largest as nu -> 0.\n")
print(f"{'nu':>6} {'alpha_num predicted':>21} {'measured peak after 1 rev':>27}")
for nu in (0.25, 0.5, 0.75, 0.9, 1.0):
    nsteps = int(round(1.0 / (nu / nx)))
    u = advance(gauss, nu, nsteps, "upwind", nx)
    a_num = C * (1.0 / nx) * (1 - nu) / 2
    print(f"{nu:>6.2f} {a_num:>21.6e} {u.max():>27.6f}")
print("\nThe peak falls monotonically as alpha_num rises. At nu = 1 the update reduces")
print("to u_j^{n+1} = u_{j-1}^n -- a pure one-cell shift -- so after nx steps the array")
print("is bit-for-bit the initial one. The peak reads 0.9983 rather than 1.0000 only")
print("because no grid point sits exactly at the Gaussian's centre.")

print()
print("=" * 78)
print("3.4  Implicit advection: unconditionally stable, and that is not the point")
print("=" * 78)
print(f"{'nu':>6} {'steps':>7} {'peak':>9} {'centroid':>10} {'exact centroid':>15} "
      f"{'phase lag':>11}")
for nu in (0.5, 1.0, 2.0, 5.0, 10.0):
    nsteps = max(1, int(round(1.0 / (nu / nx))))
    u = advance(gauss, nu, nsteps, "implicit upwind", nx)
    cen = np.sum(x * u) / np.sum(u)
    print(f"{nu:>6.1f} {nsteps:>7} {u.max():>9.4f} {cen:>10.4f} {0.3:>15.4f} "
          f"{cen - 0.3:>11.4f}")
print("\nNothing diverges at any time step. The pulse is simply annihilated -- at nu = 10")
print("the peak is 18% of its initial value after one revolution, and it has drifted")
print("almost a third of its own width downstream. Compare")
print("with the diffusion problem, where a large implicit step cost accuracy but kept")
print("the answer recognisable. Here stability buys nothing worth having.")

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))

ax[0].plot(x, gauss, "k", lw=1.3, label="exact / initial")
for scheme, c, ls in (("upwind", ORA, "-"), ("Lax-Wendroff", BLU, "-"),
                      ("Lax-Friedrichs", PUR, "--")):
    ax[0].plot(x, res31[(scheme, 0.5)], ls, color=c, lw=1.3, label=scheme)
ax[0].set_title("Gaussian after one revolution, $\\nu = 0.5$", fontsize=9)
ax[0].set_xlabel("$x$"); ax[0].set_ylabel("$u$"); ax[0].legend(fontsize=8, frameon=False)
ax[0].set_xlim(0.1, 0.5)

ax[1].plot(x, square, "k", lw=1.3, label="exact")
ax[1].plot(x, res32["upwind"], "-", color=ORA, lw=1.3, label="upwind: dissipation")
ax[1].plot(x, res32["Lax-Wendroff"], "-", color=BLU, lw=1.3, label="Lax–Wendroff: dispersion")
ax[1].axhline(0, color=GRY, lw=0.7)
ax[1].set_title("Square pulse after one revolution, $\\nu = 0.5$", fontsize=9)
ax[1].set_xlabel("$x$"); ax[1].legend(fontsize=8, frameon=False); ax[1].set_xlim(0.1, 0.55)

ax[2].plot(x, gauss, "k", lw=1.3, label="exact")
for nu, c in ((1.0, GRN), (2.0, ORA), (5.0, RED), (10.0, PUR)):
    nsteps = max(1, int(round(1.0 / (nu / nx))))
    ax[2].plot(x, advance(gauss, nu, nsteps, "implicit upwind", nx), "-",
               color=c, lw=1.3, label=f"implicit, $\\nu$ = {nu:g}")
ax[2].set_title("Implicit upwind: stable at every $\\nu$,\nand the wave is gone", fontsize=9)
ax[2].set_xlabel("$x$"); ax[2].legend(fontsize=8, frameon=False)

fig.tight_layout(); fig.savefig(f"{FIG}/time_fig3_advection.png", bbox_inches="tight")
plt.close(fig)

# --- dissipation / dispersion spectra ---------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.6))
theta = np.linspace(1e-6, np.pi, 400)              # k dx
nu = 0.5
G = {
    "upwind": 1 - nu * (1 - np.exp(-1j * theta)),
    "Lax-Wendroff": 1 - 1j * nu * np.sin(theta) - nu**2 * (1 - np.cos(theta)),
    "Lax-Friedrichs": np.cos(theta) - 1j * nu * np.sin(theta),
    "implicit upwind": 1 / (1 + nu * (1 - np.exp(-1j * theta))),
}
for (name, g), c in zip(G.items(), (ORA, BLU, PUR, RED)):
    ax[0].plot(theta / np.pi, np.abs(g), color=c, lw=1.4, label=name)
    ax[1].plot(theta / np.pi, -np.angle(g) / (nu * theta), color=c, lw=1.4, label=name)
ax[0].axhline(1, color="k", lw=1.0, ls=":")
ax[1].axhline(1, color="k", lw=1.0, ls=":")
ax[0].set_ylabel("$|g|$  — amplitude per step"); ax[0].set_title("Dissipation", fontsize=10)
ax[1].set_ylabel("$c_{num}/c$  — phase speed"); ax[1].set_title("Dispersion", fontsize=10)
for a in ax:
    a.set_xlabel(r"$k\Delta x/\pi$   (1 = the sawtooth mode)")
    a.legend(fontsize=8, frameon=False); a.grid(alpha=0.25)
ax[1].set_ylim(0, 1.35)
fig.suptitle(r"Amplification factor at $\nu = 0.5$. The exact scheme would be "
             r"$|g| = 1$ and $c_{num}/c = 1$ everywhere.", y=1.03, fontsize=9.5)
fig.tight_layout(); fig.savefig(f"{FIG}/time_fig4_spectra.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigures -> {os.path.normpath(FIG)}")

# =====================================================================
# 3.1  The modes of the advection operator
# =====================================================================
print()
print("=" * 78)
print("3.1  Where the modes are")
print("=" * 78)
NXM = 16
dxm = 1.0 / NXM
phim = 2 * np.pi * np.arange(NXM) / NXM
lam_c = -1j * (C / dxm) * np.sin(phim)
lam_u = -(C / dxm) * (1 - np.exp(-1j * phim))
print(f"nx = {NXM}, c = {C}, dx = {dxm:g}   ->   c/dx = {C/dxm:.1f}\n")
print(f"{'m':>4} {'phi':>8} {'central lambda':>22} {'upwind lambda':>24}")
for m in (0, 1, 2, 4, NXM // 2, NXM - 1):
    print(f"{m:>4} {phim[m]:>8.4f} {str(np.round(lam_c[m],2)):>22} "
          f"{str(np.round(lam_u[m],2)):>24}")
print(f"\ncentral: |lambda|_max = {np.max(np.abs(lam_c)):.1f}, all on the imaginary axis")
print(f"upwind : |lambda|_max = {np.max(np.abs(lam_u)):.1f}, max Re = "
      f"{np.max(lam_u.real):.2f}, min Re = {np.min(lam_u.real):.2f}")
print(f"\nNote the sawtooth, m = {NXM//2} (phi = pi): central gives lambda = {lam_c[NXM//2]:.1f}")
print("   -- EXACTLY ZERO. The central difference of an alternating +-1 pattern is")
print("   u_{j+1} - u_{j-1} = 0, so the scheme cannot see that mode at all: it neither")
print("   advects nor damps it. Upwind gives it lambda = %.1f, the most damped of all."
      % lam_u[NXM//2].real)
print("\nBoth scale as c/dx -- ONE power of dx, not two. Halving dx doubles |lambda|_max,")
print("and it also halves the time a wave needs to cross a cell. The two track each")
print("other, so no gap opens up: advection does not become stiff under refinement.")

fig, axm = plt.subplots(1, 2, figsize=(11.5, 3.9))
xg = (np.arange(NXM) + 0.5) * dxm
for m, c_, lbl in ((1, BLU, "$m=1$  smoothest"), (2, GRN, "$m=2$"),
                   (NXM // 2, RED, f"$m={NXM//2}$  sawtooth")):
    axm[0].plot(xg, np.cos(phim[m] * np.arange(NXM)), "-o", ms=4, color=c_, lw=1.4,
                label=f"{lbl},  $|\\lambda|$ = {abs(lam_c[m]):.1f}")
axm[0].axhline(0, color=GRY, lw=0.7); axm[0].set_ylim(-1.15, 1.8)
axm[0].set_xlabel("$x$"); axm[0].set_ylabel("mode shape (real part)")
axm[0].set_title(f"Three of the {NXM} modes, $n_x$ = {NXM}", fontsize=10)
axm[0].legend(fontsize=7.5, frameon=False, loc="upper center")

axm[1].plot(lam_c.real, lam_c.imag, "o", ms=6, color=RED, label="central")
axm[1].plot(lam_u.real, lam_u.imag, "s", ms=6, color=GRN, label="upwind")
axm[1].axhline(0, color=GRY, lw=0.6); axm[1].axvline(0, color=GRY, lw=0.6)
axm[1].set_aspect("equal")
axm[1].set_xlabel(r"Re$(\lambda)$"); axm[1].set_ylabel(r"Im$(\lambda)$")
axm[1].set_title("The two spectra: same PDE, same mesh", fontsize=10)
axm[1].set_xlim(-42, 26); axm[1].set_ylim(-30, 30)
axm[1].legend(fontsize=8, frameon=False, loc="lower left")
axm[1].annotate("central: ON the axis,\nno damping at all", xy=(0.5, 16), xytext=(4, 26),
                fontsize=7, color=RED, va="top",
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))
axm[1].annotate("the sawtooth: central gives it\n$\\lambda=0$ — it never moves",
                xy=(0.5, 0), xytext=(4, -12), fontsize=7, color=RED, va="top",
                arrowprops=dict(arrowstyle="->", color=RED, lw=0.7))
axm[1].annotate("upwind: pulled into\nthe left half-plane", xy=(-30, 11.3),
                xytext=(-41, 27), fontsize=7, color=GRN, va="top",
                arrowprops=dict(arrowstyle="->", color=GRN, lw=0.7))
fig.suptitle("Advection: the spatial scheme decides where the spectrum sits", y=1.04,
             fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/time_fig11_advection_modes.png", bbox_inches="tight")
plt.close(fig)

# =====================================================================
# 3.3  Which reference? For advection the two disagree completely
# =====================================================================
print()
print("=" * 78)
print("3.3  Spatial and temporal error do NOT separate here")
print("=" * 78)
kk = 2 * np.pi * np.fft.fftfreq(nx, d=1.0 / nx)
sym_up = -(C * nx) * (1 - np.exp(-1j * kk / nx))


def semidiscrete(t):
    return np.real(np.fft.ifft(np.fft.fft(gauss) * np.exp(sym_up * t)))


print(f"{'nu':>6} {'steps':>7} {'vs the PDE solution':>21} {'vs semi-discrete':>19}")
for nu in (0.25, 0.5, 0.9, 1.0):
    ns = int(round(1.0 / (nu / nx)))
    u = advance(gauss, nu, ns, "upwind", nx)
    print(f"{nu:>6.2f} {ns:>7} {np.max(np.abs(u - gauss)):>21.4e} "
          f"{np.max(np.abs(u - semidiscrete(1.0))):>19.4e}")
print("\nAt nu = 1 the scheme is EXACT for the PDE and maximally wrong against the")
print("semi-discrete solution. The upwind spatial operator is itself diffusive, and")
print("the forward Euler time error at nu = 1 exactly cancels that diffusion.")
print("Section 2 could isolate the temporal error; here it cannot be isolated, and")
print("the PDE solution is the reference that matters.")

# =====================================================================
# 3.6  Recommended step, across grids
# =====================================================================
print()
print("=" * 78)
print("3.6  Recommended step across grids (one revolution, T = 1, rho = 4.5)")
print("=" * 78)
RHO = 4.5
print(f"{'nx':>5} | {'scheme':<22} {'nu':>6} {'steps':>7} {'work':>9} {'max error':>11}")
print("-" * 74)
for nxr in (100, 200, 400, 800):
    xr = (np.arange(nxr) + 0.5) / nxr
    g0 = np.exp(-((xr - 0.3) / 0.06) ** 2)
    for label, scheme, nu in (("upwind, forward Euler", "upwind", 0.9),
                              ("Lax-Wendroff", "Lax-Wendroff", 0.9),
                              ("upwind, implicit", "implicit upwind", 0.9),
                              ("upwind, implicit", "implicit upwind", 5.0)):
        ns = max(1, int(round(1.0 / (nu / nxr))))
        u = advance(g0, nu, ns, scheme, nxr)
        w = ns * (RHO if "implicit" in scheme else 1.0)
        err = np.max(np.abs(u - g0))
        print(f"{nxr if label=='upwind, forward Euler' and nu==0.9 else '':>5} | "
              f"{label:<22} {nu:>6.1f} {ns:>7} {w:>9.0f} {err:>11.4f}")
    print("-" * 74)
print("Explicit at nu = 0.9 and implicit at nu = 0.9 take the SAME step; the implicit")
print("one merely costs rho times more and is slightly worse. Implicit at nu = 5 takes")
print("5x fewer steps and destroys the pulse. There is no mesh on which implicit wins.")
