"""
TIME_INTEGRATION.md sec.5 -- mixed stiffness.
   u_t + c u_x = alpha u_xx + K u (1 - u)      on a periodic [0,1]

Advection is not stiff, diffusion is stiff on a fine mesh, and the reaction is
nonlinear. Nothing here wants a single treatment.
"""
import os, time
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

NX, C, ALPHA, K, T = 400, 1.0, 0.02, 10.0, 0.2
DX = 1.0 / NX
x = (np.arange(NX) + 0.5) * DX
u0 = 0.4 * np.exp(-((x - 0.3) / 0.05) ** 2) + 0.02

dt_adv = DX / C
dt_dif = DX**2 / (2 * ALPHA)
dt_rea = 2.0 / K
print("=" * 78)
print("5.1  Three terms, three stability limits")
print("=" * 78)
print(f"nx = {NX}, dx = {DX:g}, c = {C}, alpha = {ALPHA}, K = {K}, T = {T}\n")
print(f"{'term':<26} {'explicit limit on dt':>22} {'steps to reach T':>18}")
for n, d in (("advection  dx/c", dt_adv), ("diffusion  dx^2/2alpha", dt_dif),
             ("reaction   2/K", dt_rea)):
    print(f"{n:<26} {d:>22.3e} {int(np.ceil(T/d)):>18,}")
print(f"\nThe binding limit is diffusion, by a factor of {dt_adv/dt_dif:.0f} over advection.")
print("Treating only that one term implicitly removes the whole penalty.")

# --- operators ---------------------------------------------------------------
def adv(u):                                   # first-order upwind, c > 0
    return -C * (u - np.roll(u, 1)) / DX


def dif(u):
    return ALPHA * (np.roll(u, -1) - 2 * u + np.roll(u, 1)) / DX**2


def rea(u):
    return K * u * (1 - u)


e = np.ones(NX)
L_dif = sp.diags([e, -2 * e, e], [-1, 0, 1], shape=(NX, NX), format="lil")
L_dif[0, -1] = 1.0; L_dif[-1, 0] = 1.0
L_dif = (ALPHA / DX**2) * sp.csc_matrix(L_dif)
L_adv = sp.diags([-e, e], [0, -1], shape=(NX, NX), format="lil")
L_adv[0, -1] = 1.0
L_adv = (-C / DX) * sp.csc_matrix(-L_adv)     # matches adv() above
I = sp.identity(NX, format="csc")


def run_explicit(dt):
    u, n = u0.copy(), int(np.ceil(T / dt)); dt = T / n
    for _ in range(n):
        u = u + dt * (adv(u) + dif(u) + rea(u))
    return u, n


def run_implicit(dt, tol=1e-10):
    """Backward Euler on everything; Newton on the nonlinear reaction."""
    u, n = u0.copy(), int(np.ceil(T / dt)); dt = T / n
    A_lin = I - dt * (L_dif + L_adv)
    for _ in range(n):
        un, w = u.copy(), u.copy()
        for _ in range(30):                   # Newton
            F = w - un - dt * (L_dif @ w + L_adv @ w + rea(w))
            if np.max(np.abs(F)) < tol:
                break
            J = A_lin - dt * sp.diags(K * (1 - 2 * w))
            w = w - spla.spsolve(J.tocsc(), F)
        u = w
    return u, n


def run_imex(dt):
    """Diffusion implicit, advection and reaction explicit."""
    u, n = u0.copy(), int(np.ceil(T / dt)); dt = T / n
    lu = spla.splu((I - dt * L_dif).tocsc())
    for _ in range(n):
        u = lu.solve(u + dt * (adv(u) + rea(u)))
    return u, n


# reference: explicit, far below every limit
ref, _ = run_explicit(0.02 * min(dt_adv, dt_dif, dt_rea))

print()
print("=" * 78)
print("5.2  Cost to reach T, at each scheme's own largest usable step")
print("=" * 78)
print(f"{'scheme':<34} {'dt':>10} {'steps':>9} {'wall':>9} {'max error':>12}")
rows = []
for label, fn, dt in (("fully explicit (forward Euler)", run_explicit, 0.9 * dt_dif),
                      ("fully implicit (backward Euler)", run_implicit, 0.9 * dt_adv),
                      ("IMEX: diffusion implicit", run_imex, 0.9 * dt_adv)):
    t0 = time.perf_counter(); u, n = fn(dt); wall = time.perf_counter() - t0
    err = np.max(np.abs(u - ref))
    print(f"{label:<34} {T/n:>10.2e} {n:>9,} {wall:>8.3f}s {err:>12.3e}")
    rows.append((label, T / n, n, wall, err, u))
print("\nIMEX takes the same step as the fully implicit scheme and pays for one banded")
print("solve instead of a Newton iteration; the fully explicit scheme pays for neither")
print("but takes 16 times as many steps.")

print()
print("=" * 78)
print("5.2b  The same comparison at EQUAL ACCURACY, which is the fair one")
print("=" * 78)
TARGET = 1e-3
print(f"largest dt for which max error <= {TARGET:g}, found by halving:\n")
print(f"{'scheme':<34} {'dt':>10} {'steps':>9} {'wall':>9} {'error':>12}")
fair = []
for label, fn, dt0 in (("fully explicit (forward Euler)", run_explicit, 0.9 * dt_dif),
                       ("fully implicit (backward Euler)", run_implicit, 0.9 * dt_adv),
                       ("IMEX: diffusion implicit", run_imex, 0.9 * dt_adv)):
    dt = dt0
    while True:
        u, n = fn(dt)
        if np.max(np.abs(u - ref)) <= TARGET or dt < 1e-6:
            break
        dt /= 2
    t0 = time.perf_counter(); u, n = fn(dt); wall = time.perf_counter() - t0
    err = np.max(np.abs(u - ref))
    print(f"{label:<34} {T/n:>10.2e} {n:>9,} {wall:>8.3f}s {err:>12.3e}")
    fair.append((label, T / n, n, wall, err))
print("\nThe stability advantage does not survive. Forward Euler, IMEX-Euler and")
print("backward Euler are ALL first order, so hitting a fixed error target forces")
print("them to nearly the same step, and the explicit scheme -- which pays for no")
print("solve at all -- is then competitive or better.")
print("\nThe lesson is not that IMEX is useless. It is that a stability-limited")
print("comparison flatters implicit methods, and that the real case for IMEX needs")
print("either a higher-order IMEX pair, or a problem where the stiff mode is")
print("genuinely uninteresting so that its error does not matter.")

# =====================================================================
# 5.3  Splitting order, with both sub-solves exact
# =====================================================================
print()
print("=" * 78)
print("5.3  Godunov against Strang splitting, with exact sub-solves")
print("=" * 78)
print("The linear part (advection + diffusion) is solved exactly in Fourier space and")
print("the logistic reaction exactly in closed form, so the ONLY error left is the")
print("splitting error itself.\n")

kk = 2 * np.pi * np.fft.fftfreq(NX, d=DX)
sym_lin = (-C / DX) * (1 - np.exp(-1j * kk * DX)) \
          + (ALPHA / DX**2) * (np.exp(1j * kk * DX) - 2 + np.exp(-1j * kk * DX))


def exact_lin(u, dt):
    return np.real(np.fft.ifft(np.fft.fft(u) * np.exp(sym_lin * dt)))


def exact_rea(u, dt):                        # logistic ODE, closed form
    E = np.exp(K * dt)
    return u * E / (1 - u + u * E)


def split(dt, kind):
    u, n = u0.copy(), int(np.ceil(T / dt)); dt = T / n
    for _ in range(n):
        if kind == "godunov":
            u = exact_rea(exact_lin(u, dt), dt)
        else:                                 # Strang
            u = exact_lin(u, dt / 2)
            u = exact_rea(u, dt)
            u = exact_lin(u, dt / 2)
    return u


ref_split = split(1e-6, "strang")             # unsplit limit
print(f"{'dt':>10} {'Godunov error':>16} {'p':>6} {'Strang error':>16} {'p':>6}")
prev = {"godunov": None, "strang": None}
dts = [4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4]
conv = {"godunov": [], "strang": []}
for dt in dts:
    row = f"{dt:>10.2e} "
    for kind in ("godunov", "strang"):
        err = np.max(np.abs(split(dt, kind) - ref_split))
        conv[kind].append(err)
        p = "" if prev[kind] is None else f"{np.log(prev[kind]/err)/np.log(2):.2f}"
        row += f"{err:>16.3e} {p:>6}"
        prev[kind] = err
    print(row)
print("\nGodunov first order, Strang second, for the price of splitting the cheap")
print("half-step in two. Splitting is not free accuracy -- but Strang makes it cheap.")

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))

ax[0].plot(x, u0, "k--", lw=1.1, label="$t = 0$")
ax[0].plot(x, ref, "k", lw=1.6, label="reference")
for (label, dt, n, wall, err, u), c in zip(rows, (RED, BLU, GRN)):
    ax[0].plot(x, u, "-", color=c, lw=1.2, alpha=0.9, label=label.split(" (")[0])
ax[0].set_xlabel("$x$"); ax[0].set_ylabel("$u$")
ax[0].set_title(f"Solution at $T$ = {T}", fontsize=9)
ax[0].legend(fontsize=7, frameon=False)

labs = [r[0].split(" (")[0] for r in rows]
ax[1].bar(range(3), [r[2] for r in rows], color=[RED, BLU, GRN], alpha=0.85)
ax[1].set_yscale("log"); ax[1].set_xticks(range(3))
ax[1].set_xticklabels(["explicit", "implicit", "IMEX"], fontsize=8)
ax[1].set_ylabel("steps to reach $T$")
for i, r in enumerate(rows):
    ax[1].text(i, r[2] * 1.3, f"{r[2]:,}\n{r[3]:.2f}s", ha="center", fontsize=7.5)
ax[1].set_title("Step count and wall time", fontsize=9)
ax[1].set_ylim(top=max(r[2] for r in rows) * 6)

ax[2].loglog(dts, conv["godunov"], "o-", color=RED, label="Godunov (1st order)")
ax[2].loglog(dts, conv["strang"], "s-", color=GRN, label="Strang (2nd order)")
d = np.array(dts)
ax[2].loglog(d, conv["godunov"][0] * d / d[0], ":", color=GRY, lw=1)
ax[2].loglog(d, conv["strang"][0] * (d / d[0])**2, ":", color=GRY, lw=1)
ax[2].set_xlabel(r"$\Delta t$"); ax[2].set_ylabel("splitting error")
ax[2].set_title("Splitting error, exact sub-solves", fontsize=9)
ax[2].legend(fontsize=8, frameon=False); ax[2].grid(alpha=0.25, which="both")

fig.tight_layout(); fig.savefig(f"{FIG}/time_fig6_imex.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigures -> {os.path.normpath(FIG)}")

# =====================================================================
# 5.0 / 5.1  The three mode clusters
# =====================================================================
print()
print("=" * 78)
print("5.1b  The spectrum is a SUM of three, with different directions and scalings")
print("=" * 78)
phi = 2 * np.pi * np.arange(NX) / NX
lam_adv = -(C / DX) * (1 - np.exp(-1j * phi))          # upwind
lam_dif = (ALPHA / DX**2) * (2 * np.cos(phi) - 2)      # central
lam_rea_0 = K                                          # R'(u) at u = 0  (GROWTH)
lam_rea_1 = -K                                         # R'(u) at u = 1  (decay)
lam_tot = lam_adv + lam_dif
print(f"nx = {NX}, dx = {DX:g}, c = {C}, alpha = {ALPHA}, K = {K}\n")
print(f"{'branch':<24} {'|lambda|_max':>14} {'direction':>26} {'scaling':>12}")
print(f"{'advection (upwind)':<24} {np.max(np.abs(lam_adv)):>14.1f} "
      f"{'circle, Re < 0':>26} {'1/dx':>12}")
print(f"{'diffusion':<24} {np.max(np.abs(lam_dif)):>14.1f} "
      f"{'negative real axis':>26} {'1/dx^2':>12}")
print(f"{'reaction, linearised':<24} {K:>14.1f} "
      f"{'real: +K at u=0, -K at u=1':>26} {'1 (no dx)':>12}")
print(f"\ncombined advection+diffusion: |lambda|_max = {np.max(np.abs(lam_tot)):.1f}, "
      f"max Re = {np.max(lam_tot.real):.1f}, min Re = {np.min(lam_tot.real):.1f}")
print("\nNote the reaction eigenvalue is POSITIVE where u is small: R(u)=Ku(1-u) has")
print(f"R'(0) = +{K:g}, so the true solution GROWS there. This is the one place in the")
print("document where |g| <= 1 is the wrong stability requirement (sec.0.5).")

fig, ax5 = plt.subplots(1, 2, figsize=(11.5, 4.0))
ax5[0].plot(lam_adv.real, lam_adv.imag, "o", ms=5, color=BLU, label="advection alone")
ax5[0].plot(lam_dif.real, lam_dif.imag, "s", ms=5, color=RED, label="diffusion alone")
ax5[0].plot([lam_rea_0, lam_rea_1], [0, 0], "^", ms=9, color=GRN,
            label="reaction ($\\pm K$)")
ax5[0].axhline(0, color=GRY, lw=0.6); ax5[0].axvline(0, color=GRY, lw=0.6)
ax5[0].set_xlabel(r"Re$(\lambda)$"); ax5[0].set_ylabel(r"Im$(\lambda)$")
ax5[0].set_title("The three branches separately", fontsize=10)
ax5[0].legend(fontsize=7.5, frameon=False, loc="upper left")
ax5[0].set_yscale("symlog", linthresh=100); ax5[0].set_xscale("symlog", linthresh=100)

ax5[1].plot(lam_tot.real, lam_tot.imag, "o", ms=5, color=PUR, label="advection + diffusion")
ax5[1].plot([lam_rea_0], [0], "^", ms=9, color=GRN, label="reaction, $+K$ (growth)")
ax5[1].axhline(0, color=GRY, lw=0.6); ax5[1].axvline(0, color=GRY, lw=0.6)
ax5[1].set_xlabel(r"Re$(\lambda)$"); ax5[1].set_ylabel(r"Im$(\lambda)$")
ax5[1].set_title("Combined — diffusion dominates the extent", fontsize=10)
ax5[1].legend(fontsize=7.5, frameon=False, loc="upper left")
ax5[1].set_xscale("symlog", linthresh=100)
fig.suptitle("Three physical terms, three clusters: no single treatment suits all of them",
             y=1.03, fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/time_fig13_imex_modes.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigure -> {os.path.normpath(FIG)}/time_fig13_imex_modes.png")
