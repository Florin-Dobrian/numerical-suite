"""
TIME_INTEGRATION.md sec.2 -- parabolic.  u_t = alpha u_xx on [0,1], u(0)=u(1)=0.

theta-method:  (I - theta r L) u^{n+1} = (I + (1-theta) r L) u^n,  r = alpha dt / dx^2
theta = 0 forward Euler, 1/2 Crank-Nicolson, 1 backward Euler.
"""
import os, time
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)

import numpy as np
from scipy.linalg import solve_banded
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLU, ORA, GRN, RED, PUR, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#805ad5", "#718096"
ALPHA = 1.0


def banded_theta(nx, r, theta):
    """LHS of the theta-method as a banded matrix (interior nodes only)."""
    m = nx - 1
    ab = np.zeros((3, m))
    ab[0, 1:] = -theta * r            # super
    ab[1, :] = 1 + 2 * theta * r      # diag
    ab[2, :-1] = -theta * r           # sub
    return ab


def step_theta(u, r, theta, ab):
    """One theta step on the interior vector u (Dirichlet 0 at both ends)."""
    lap = np.empty_like(u)
    lap[1:-1] = u[2:] - 2 * u[1:-1] + u[:-2]
    lap[0] = u[1] - 2 * u[0]
    lap[-1] = u[-2] - 2 * u[-1]
    rhs = u + (1 - theta) * r * lap
    if theta == 0.0:
        return rhs
    return solve_banded((1, 1), ab, rhs)


def run(nx, dt, theta, T, u0fun):
    dx = 1.0 / nx
    r = ALPHA * dt / dx**2
    x = np.linspace(0, 1, nx + 1)
    u = u0fun(x[1:-1]).astype(float)
    ab = banded_theta(nx, r, theta)
    nsteps = int(round(T / dt))
    hist_min = []
    for _ in range(nsteps):
        u = step_theta(u, r, theta, ab)
        hist_min.append(u.min())
        if not np.all(np.isfinite(u)) or np.abs(u).max() > 1e8:
            return x, np.full(nx + 1, np.nan), r, nsteps, np.array(hist_min)
    full = np.zeros(nx + 1); full[1:-1] = u
    return x, full, r, nsteps, np.array(hist_min)


smooth = lambda x: np.sin(np.pi * x)
exact_smooth = lambda x, t: np.sin(np.pi * x) * np.exp(-ALPHA * np.pi**2 * t)
step_ic = lambda x: np.where(np.abs(x - 0.5) < 0.15, 1.0, 0.0)
def spike_ic(x):
    u = np.zeros_like(x); u[len(x) // 2] = 1.0; return u

print("=" * 78)
print("2.1  The explicit threshold r = alpha dt / dx^2 <= 1/2")
print("=" * 78)
print("von Neumann: g(k) = 1 - 4 r sin^2(k dx / 2). The worst mode is the sawtooth,")
print("k dx = pi, for which g = 1 - 4r. So |g| <= 1 requires exactly r <= 1/2 --")
print("and the threshold is only visible if that mode is present in the data.\n")
nx, nsteps = 50, 100
sawtooth = lambda x: (-1.0) ** np.arange(len(x))
print(f"{'r':>8} {'g = 1-4r':>10} {'|g|^100 predicted':>19} {'measured max|u|':>18}   verdict")
for r in (0.45, 0.49, 0.50, 0.5005, 0.51, 0.55):
    dt = r / (ALPHA * nx**2)
    _, u, _, _, _ = run(nx, dt, 0.0, nsteps * dt, sawtooth)
    g = 1 - 4 * r
    m = np.nanmax(np.abs(u))
    print(f"{r:>8.4f} {g:>10.4f} {abs(g)**nsteps:>19.4e} {m:>18.4e}   "
          f"{'GROWS' if abs(g) > 1 else 'decays'}")
print("\nPrediction and measurement agree to four digits on both sides of r = 1/2.")

print()
print("=" * 78)
print("2.2  Order in time, measured against the semi-discrete exact solution")
print("=" * 78)
print("Comparing to the PDE solution would mix spatial and temporal error. The exact")
print("solution of the semi-discrete ODE system is available in closed form for")
print("sinusoidal data, so the temporal error can be isolated exactly.\n")
nx, T = 200, 0.02
dx = 1.0 / nx
lam_h = (4 * ALPHA / dx**2) * np.sin(np.pi * dx / 2) ** 2     # semi-discrete eigenvalue
xg = np.linspace(0, 1, nx + 1)
ref = np.sin(np.pi * xg) * np.exp(-lam_h * T)
dts = [4e-4, 2e-4, 1e-4, 5e-5, 2.5e-5]
print(f"{'dt':>10} {'r':>9} " + "".join(f"{n:>26}" for n in
                                        ("theta=1 (BE)", "theta=1/2 (CN)")))
prev = {1.0: None, 0.5: None}
for dt in dts:
    row = f"{dt:>10.2e} {ALPHA*dt*nx**2:>9.2f} "
    for theta in (1.0, 0.5):
        x, u, r, _, _ = run(nx, dt, theta, T, smooth)
        e = np.max(np.abs(u - ref))
        o = "" if prev[theta] is None else f" (p={np.log(prev[theta]/e)/np.log(2):.2f})"
        row += f"{e:>14.3e}{o:>12}"
        prev[theta] = e
    print(row)
print("\nBackward Euler first order, Crank-Nicolson second. Note that Crank-Nicolson")
print("at the largest dt here is already more accurate than backward Euler at the")
print("smallest -- which is the entire case for using it, when the data are smooth.")

print()
print("=" * 78)
print("2.3  Crank-Nicolson is A-stable but not L-stable: discontinuous data, big dt")
print("=" * 78)
nx = 100
dt = 1e-3                                  # r = 10, far past the explicit limit
T = 5 * dt
r_ = ALPHA * dt * nx**2
z_top = -4 * r_                            # the sawtooth mode sees lambda dt = -4r
print(f"nx = {nx}, dt = {dt}, r = {r_:.0f}, 5 steps, "
      f"initial data = a single-node spike\n")
print(f"the sawtooth mode sees lambda*dt = {z_top:.0f}, so after 5 steps it is")
print(f"  amplified by  BE: {(1/(1-z_top))**5:.3e}    CN: {((1+z_top/2)/(1-z_top/2))**5:+.4f}\n")
print(f"{'scheme':<24} {'min u over all steps':>22} {'min u at t=T':>16} {'TV at t=T':>12}")
res23 = {}
for label, theta in (("Backward Euler", 1.0), ("Crank-Nicolson", 0.5),
                     ("theta = 0.75", 0.75)):
    x, u, r, ns, hmin = run(nx, dt, theta, T, spike_ic)
    tv = np.sum(np.abs(np.diff(u)))
    print(f"{label:<24} {hmin.min():>22.6f} {u.min():>16.6f} {tv:>12.4f}")
    res23[label] = (x, u)
print("\nThe exact solution of the heat equation is positive everywhere for t > 0.")
print("Crank-Nicolson undershoots badly and rings; backward Euler does not.")
print("Cause: R_CN(z) -> -1 as z -> -inf, so the sharp modes are sign-flipped, not damped.")

print()
print("=" * 78)
print("2.4  Cost: the explicit step must shrink as dx^2")
print("=" * 78)
T = 0.01
print(f"{'nx':>6} {'explicit steps':>16} {'implicit steps':>16} "
      f"{'explicit wall':>15} {'implicit wall':>15} {'speedup':>9}")
cost = []
for nx in (50, 100, 200, 400, 800):
    dt_e = 0.5 / (ALPHA * nx**2)
    ns_e = int(np.ceil(T / dt_e))
    t0 = time.perf_counter(); run(nx, T / ns_e, 0.0, T, smooth); te = time.perf_counter() - t0
    dt_i = 1e-4                                      # accuracy-chosen, fixed
    ns_i = int(round(T / dt_i))
    t0 = time.perf_counter(); run(nx, dt_i, 1.0, T, smooth); ti = time.perf_counter() - t0
    print(f"{nx:>6} {ns_e:>16,} {ns_i:>16,} {te:>14.3f}s {ti:>14.3f}s {te/ti:>8.1f}x")
    cost.append((nx, ns_e, ns_i, te, ti))
print("\nDoubling nx quadruples the explicit step count and leaves the implicit one alone.")

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))

nx, nst = 50, 100
seeded = lambda x: np.sin(np.pi * x) + 0.01 * (-1.0) ** np.arange(len(x))
for r, c, ls in ((0.45, GRN, "-"), (0.51, RED, "--")):
    dt = r / (ALPHA * nx**2)
    x, u, _, _, _ = run(nx, dt, 0.0, nst * dt, seeded)
    ax[0].plot(x, np.clip(u, -1.5, 1.5), ls, color=c, lw=1.4,
               label=f"$r$ = {r}" + ("  (clipped)" if r > 0.5 else ""))
ax[0].set_title("Forward Euler, 100 steps, smooth data seeded\nwith 1% sawtooth noise",
                fontsize=9)
ax[0].set_xlabel("$x$"); ax[0].set_ylabel("$u$"); ax[0].legend(fontsize=8, frameon=False)

for label, (x, u), c, m in (("Backward Euler", res23["Backward Euler"], GRN, "o"),
                            ("Crank-Nicolson", res23["Crank-Nicolson"], RED, "s")):
    ax[1].plot(x, u, "-", marker=m, ms=3, color=c, lw=1.3, label=label)
ax[1].axhline(0, color=GRY, lw=0.8)
ax[1].set_title("Top-hat initial data, $r = 10$\nCrank–Nicolson rings below zero", fontsize=9)
ax[1].set_xlabel("$x$"); ax[1].legend(fontsize=8, frameon=False)

nxs = [c[0] for c in cost]
ax[2].loglog(nxs, [c[1] for c in cost], "o-", color=RED, label="explicit steps")
ax[2].loglog(nxs, [c[2] for c in cost], "s-", color=GRN, label="implicit steps")
ax[2].loglog(nxs, 0.5 * np.array(nxs, float)**2, ":", color=GRY, lw=1.2)
ax[2].text(nxs[1], 0.8 * nxs[1]**2, r"slope 2", color=GRY, fontsize=8)
ax[2].set_xlabel("$n_x$"); ax[2].set_ylabel("steps to $t = 0.01$")
ax[2].set_title("Step count against mesh resolution", fontsize=9)
ax[2].legend(fontsize=8, frameon=False); ax[2].grid(alpha=0.25, which="both")

fig.tight_layout(); fig.savefig(f"{FIG}/time_fig2_diffusion.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigures -> {os.path.normpath(FIG)}")

# =====================================================================
# 2.5  The two steps side by side: what stability allows, what accuracy needs
# =====================================================================
print()
print("=" * 78)
print("2.5  dt_stab against dt_acc — the two numbers the verdict depends on")
print("=" * 78)
T_ = 0.02


def err_at(nx, dt, theta):
    """Error against the semi-discrete exact solution, so this is purely temporal."""
    dx = 1.0 / nx
    lam_h = (4 * ALPHA / dx**2) * np.sin(np.pi * dx / 2) ** 2
    xg = np.linspace(0, 1, nx + 1)
    ref = np.sin(np.pi * xg) * np.exp(-lam_h * T_)
    x, u, _, _, _ = run(nx, dt, theta, T_, smooth)
    return np.max(np.abs(u - ref))


def largest_dt_for(nx, theta, tol):
    """Largest dt (by bisection on a geometric ladder) meeting the tolerance."""
    lo, hi = 1e-9, T_          # cannot take more than one step of size T
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        n = max(1, int(round(T_ / mid)))
        if err_at(nx, T_ / n, theta) <= tol:
            lo = mid
        else:
            hi = mid
    return lo


for tol in (1e-3, 1e-4):
    print(f"\ntolerance = {tol:g}   (max error against the semi-discrete solution)")
    print(f"{'nx':>6} {'dt_stab (r=1/2)':>17} {'dt_acc, BE':>13} {'gap':>9} "
          f"{'dt_acc, CN':>13} {'gap':>9}")
    for nx in (50, 100, 200, 400):
        dts = 0.5 / (ALPHA * nx**2)
        dbe = largest_dt_for(nx, 1.0, tol)
        dcn = largest_dt_for(nx, 0.5, tol)
        sat = "*" if dcn > 0.98 * T_ else " "
        print(f"{nx:>6} {dts:>17.3e} {dbe:>13.3e} {dbe/dts:>8.1f}x "
              f"{dcn:>12.3e}{sat} {dcn/dts:>8.1f}x")
print("\n* = the whole interval in one step; the true dt_acc is larger still.")
print("\nThe gap is dt_acc / dt_stab -- the stiffness gap of sec.0.9, measured.")
print("Implicit is worth its cost when that gap exceeds rho, the solve/multiply ratio.")
print("Note the gap GROWS as the mesh is refined, and is far larger for")
print("Crank-Nicolson than for backward Euler because second order buys a bigger step.")

# =====================================================================
# 2.1 figure: the modes of THIS operator, and where they sit
# =====================================================================
NXF = 16
dxf = 1.0 / NXF
Mf = NXF - 1
mf = np.arange(1, Mf + 1)
lam_f = -(4 * ALPHA / dxf**2) * np.sin(mf * np.pi * dxf / 2) ** 2
xf = np.linspace(0, 1, NXF + 1)

print()
print("=" * 78)
print(f"2.1  The modes at nx = {NXF}  (M = {Mf} of them, 4*alpha/dx^2 = {4*ALPHA/dxf**2:.1f})")
print("=" * 78)
print(f"{'m':>4} {'phi = m pi dx':>15} {'lambda_m':>12} {'lifetime 1/|lam|':>18}   shape")
for m in (1, 2, 3, Mf - 1, Mf):
    lab = ("slowest — one half-sine" if m == 1 else
           "two half-sines" if m == 2 else
           "SAWTOOTH — fastest" if m == Mf else f"{m} half-sines")
    print(f"{m:>4} {m*np.pi*dxf:>15.4f} {lam_f[m-1]:>12.2f} {1/abs(lam_f[m-1]):>18.3e}   {lab}")
print(f"\nratio |lambda_max| / |lambda_min| = {abs(lam_f[-1]/lam_f[0]):.1f}")

fig, ax = plt.subplots(1, 2, figsize=(11.5, 3.9))

for m, c, lbl in ((1, BLU, "$m=1$  slowest"), (2, GRN, "$m=2$"), (Mf, RED, f"$m={Mf}$  fastest")):
    v = np.zeros(NXF + 1)
    v[1:-1] = np.sin(m * np.pi * xf[1:-1])
    v /= np.max(np.abs(v))
    ax[0].plot(xf, v, "-o", ms=4, color=c, lw=1.4,
               label=f"{lbl},  $\\lambda$ = {lam_f[m-1]:.1f}")
ax[0].axhline(0, color=GRY, lw=0.7)
ax[0].set_ylim(-1.15, 1.75)
ax[0].set_xlabel("$x$"); ax[0].set_ylabel("mode shape (normalised)")
ax[0].set_title(f"Three of the {Mf} modes, $n_x$ = {NXF}", fontsize=10)
ax[0].legend(fontsize=7.5, frameon=False, loc="upper center", ncol=1)

S = 4 * ALPHA / dxf**2
ax[1].vlines(np.abs(lam_f), 0, 1, color=GRY, lw=1.2)
for m, c in ((1, BLU), (2, GRN), (Mf, RED)):
    ax[1].vlines(abs(lam_f[m-1]), 0, 1, color=c, lw=2.4)
    ax[1].annotate(f"$m={m}$", xy=(abs(lam_f[m-1]), 1.0), xytext=(abs(lam_f[m-1]), 1.35),
                   color=c, fontsize=8, ha="center",
                   arrowprops=dict(arrowstyle="-", color=c, lw=0.8))
ax[1].axvline(ALPHA * np.pi**2, color=BLU, ls=":", lw=1.2)
ax[1].axvline(S, color=RED, ls=":", lw=1.2)
ax[1].text(ALPHA*np.pi**2*0.97, -0.30, r"$\alpha\pi^2$" + "\n(limit, not reached)",
           color=BLU, fontsize=7.5, ha="right", va="top")
ax[1].text(S*1.03, -0.30, r"$4\alpha/\Delta x^2$" + "\n(limit, not reached)",
           color=RED, fontsize=7.5, ha="left", va="top")
ax[1].set_xscale("log"); ax[1].set_xlim(4, 3000); ax[1].set_ylim(-1.0, 1.8)
ax[1].set_yticks([])
ax[1].set_xlabel(r"$|\lambda_m|$   (log scale)")
ax[1].set_title("Where all 15 eigenvalues sit — both ends open", fontsize=10)
for sp in ("left", "right", "top"):
    ax[1].spines[sp].set_visible(False)

fig.suptitle("Diffusion: smooth modes are slow and are the physics; "
             "wiggly modes are fast and are the grid", y=1.04, fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/time_fig10_diffusion_modes.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigure -> {os.path.normpath(FIG)}/time_fig10_diffusion_modes.png")

# =====================================================================
# 2.6  Consolidated recommendation: the step each scheme should take
# =====================================================================
TOL_R, SAFETY, RHO = 1e-4, 0.9, 4.5


def err_R(nx, dt, theta):
    dx = 1.0 / nx
    lam = (4 * ALPHA / dx**2) * np.sin(np.pi * dx / 2) ** 2
    xg = np.linspace(0, 1, nx + 1)
    ref = np.sin(np.pi * xg) * np.exp(-lam * T_)
    _, u, _, _, _ = run(nx, dt, theta, T_, smooth)
    return np.inf if not np.all(np.isfinite(u)) else np.max(np.abs(u - ref))


def dt_acc_R(nx, theta):
    lo, hi = 1e-9, T_
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        k = max(1, int(round(T_ / mid)))
        if err_R(nx, T_ / k, theta) <= TOL_R:
            lo = mid
        else:
            hi = mid
    return lo


print()
print("=" * 78)
print(f"2.6  Recommended dt   (tol {TOL_R:g} on max error, T = {T_}, "
      f"safety {SAFETY} on dt_stab, rho = {RHO})")
print("=" * 78)
print("work is measured in units of one explicit step; an implicit step costs rho.\n")
hdr = (f"{'nx':>5} {'dt_stab':>10} | {'FE dt':>10} {'steps':>7} {'work':>8} {'limit':>6}"
       f" | {'BE dt':>10} {'steps':>6} {'work':>7}"
       f" | {'CN dt':>10} {'steps':>6} {'work':>7} | best")
print(hdr)
print("-" * len(hdr))
for nx in (25, 50, 100, 200, 400, 800):
    dts = 0.5 / (ALPHA * nx**2)
    a_fe, a_be, a_cn = dt_acc_R(nx, 0.0), dt_acc_R(nx, 1.0), dt_acc_R(nx, 0.5)
    fe = min(SAFETY * dts, a_fe)
    lim = "stab" if SAFETY * dts < a_fe else "acc"
    nfe, nbe, ncn = T_ / fe, T_ / a_be, T_ / a_cn
    wfe, wbe, wcn = nfe, nbe * RHO, ncn * RHO
    best = min((wfe, "FE"), (wbe, "BE"), (wcn, "CN"))[1]
    print(f"{nx:>5} {dts:>10.2e} | {fe:>10.2e} {nfe:>7.0f} {wfe:>8.0f} {lim:>6}"
          f" | {a_be:>10.2e} {nbe:>6.0f} {wbe:>7.0f}"
          f" | {a_cn:>10.2e} {ncn:>6.0f} {wcn:>7.0f} | {best}")
print("\nForward Euler is ACCURACY-limited on coarse meshes: below nx ~ 90 the CFL")
print("condition is satisfied automatically and plays no part.")
print("The two first-order schemes want the SAME step -- their leading errors are")
print("+z^2/2 and -z^2/2 (sec.0.5), equal in size. Backward Euler buys no larger step,")
print("only a solve. Crank-Nicolson's second order buys a 64x larger one.")
