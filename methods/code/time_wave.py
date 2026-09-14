"""
TIME_INTEGRATION.md sec.4 -- second-order hyperbolic.  u_tt = c^2 u_xx.

This is the structural-dynamics problem of COMPUTATIONAL.md sec.2.3, and the right
question is not pointwise error but whether the scheme conserves energy.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(_HERE, "..", "figs")
os.makedirs(FIG, exist_ok=True)

import numpy as np
RHO_W = 4.5
from scipy.linalg import solve_banded
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLU, ORA, GRN, RED, PUR, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#805ad5", "#718096"
Cw = 1.0

# =====================================================================
# 4.1  One oscillator: period elongation and amplitude decay, exactly
# =====================================================================
def newmark_amp(gamma, beta, Om):
    """Amplification matrix of Newmark for  u'' + w^2 u = 0,  Om = w*dt."""
    A = np.zeros((3, 3))
    for k, s in enumerate(np.eye(3)):                 # states (u, dt*v, dt^2*a)
        u, dv, da = s
        num = -Om**2 * (u + dv + (0.5 - beta) * da)
        na = num / (1 + beta * Om**2)
        nu_ = u + dv + (0.5 - beta) * da + beta * na
        ndv = dv + (1 - gamma) * da + gamma * na
        A[:, k] = [nu_, ndv, na]
    return A


def be_amp(Om):
    """Backward Euler on the first-order system (u, dt*v)."""
    return np.linalg.inv(np.array([[1.0, -1.0], [Om**2, 1.0]]))


def spectral(A):
    ev = np.linalg.eigvals(A)
    ev = ev[np.argsort(-np.abs(np.imag(ev)))]         # the oscillatory pair
    lam = ev[0]
    rho, Omb = np.abs(lam), abs(np.angle(lam))
    return rho, Omb


print("=" * 78)
print("4.1  One oscillator, u'' + w^2 u = 0. Period error and amplitude decay")
print("=" * 78)
print("Measured from the eigenvalues of each scheme's amplification matrix, so these")
print("are exact properties of the scheme, not the outcome of a particular run.\n")
SCH = {
    "central difference (explicit)": lambda Om: newmark_amp(0.5, 0.0, Om),
    "Newmark trapezoidal (1/2,1/4)": lambda Om: newmark_amp(0.5, 0.25, Om),
    "Newmark linear accel (1/2,1/6)": lambda Om: newmark_amp(0.5, 1 / 6, Om),
    "Newmark damped (0.6,0.3025)": lambda Om: newmark_amp(0.6, 0.3025, Om),
    "backward Euler": be_amp,
}
for frac in (0.01, 0.05, 0.10):
    Om = 2 * np.pi * frac                              # dt / T = frac
    print(f"dt / T = {frac:.2f}   (Omega = w*dt = {Om:.4f})")
    print(f"  {'scheme':<32} {'|lambda|':>10} {'period error':>14} {'decay/cycle':>13}")
    for name, f in SCH.items():
        rho, Omb = spectral(f(Om))
        Tn_over_T = Om / Omb
        decay = 1 - rho ** (2 * np.pi / Omb)
        print(f"  {name:<32} {rho:>10.6f} {100*(Tn_over_T-1):>13.3f}% {100*decay:>12.3f}%")
    print()
print("|lambda| = 1 exactly means no amplitude decay whatsoever: the scheme neither")
print("adds nor removes energy. Only the two gamma = 1/2 Newmark members and the")
print("explicit central difference have it. Backward Euler destroys the oscillation.")

print()
print("=" * 78)
print("4.2  The explicit stability limit for the second-order system")
print("=" * 78)
print(f"{'Omega = w dt':>13} {'|lambda| central diff':>23}   verdict")
for Om in (1.5, 1.9, 1.99, 2.0, 2.01, 2.1):
    rho = max(np.abs(np.linalg.eigvals(newmark_amp(0.5, 0.0, Om))))
    print(f"{Om:>13.3f} {rho:>23.6f}   {'UNSTABLE' if rho > 1 + 1e-9 else 'stable'}")
print("\nThe limit is Omega = w_max dt <= 2. With w_max = 2c/dx from the discrete")
print("Laplacian, that is exactly the CFL condition  c dt / dx <= 1.")

# =====================================================================
# 4.3  The PDE: energy drift over many periods
# =====================================================================
nx = 200
dx = 1.0 / nx
xg = np.linspace(0, 1, nx + 1)
m = nx - 1                                            # interior unknowns


def lap(u):
    out = np.empty_like(u)
    out[1:-1] = u[2:] - 2 * u[1:-1] + u[:-2]
    out[0] = u[1] - 2 * u[0]
    out[-1] = u[-2] - 2 * u[-1]
    return out / dx**2


def energy(u, v):
    full = np.zeros(nx + 1); full[1:-1] = u
    return 0.5 * dx * np.sum(v**2) + 0.5 * Cw**2 * np.sum(np.diff(full) ** 2) / dx


def run_newmark(u0, v0, dt, nsteps, gamma, beta, every=1):
    u, v = u0.copy(), v0.copy()
    a = Cw**2 * lap(u)
    ab = np.zeros((3, m))                             # (I + beta dt^2 c^2 (-L))
    k = beta * dt**2 * Cw**2 / dx**2
    ab[0, 1:] = -k; ab[1, :] = 1 + 2 * k; ab[2, :-1] = -k
    E = []
    for n in range(nsteps):
        up = u + dt * v + dt**2 * (0.5 - beta) * a
        rhs = Cw**2 * lap(up)
        a_new = rhs if beta == 0.0 else solve_banded((1, 1), ab, rhs)
        u = up + beta * dt**2 * a_new
        v = v + dt * ((1 - gamma) * a + gamma * a_new)
        a = a_new
        if n % every == 0:
            E.append(energy(u, v))
    return u, v, np.array(E)


def run_be(u0, v0, dt, nsteps, every=1):
    """Backward Euler on the first-order system."""
    u, v = u0.copy(), v0.copy()
    ab = np.zeros((3, m))
    k = dt**2 * Cw**2 / dx**2
    ab[0, 1:] = -k; ab[1, :] = 1 + 2 * k; ab[2, :-1] = -k
    E = []
    for n in range(nsteps):
        u = solve_banded((1, 1), ab, u + dt * v)
        v = v + dt * Cw**2 * lap(u)
        if n % every == 0:
            E.append(energy(u, v))
    return u, v, np.array(E)


print()
print("=" * 78)
print("4.3  Energy over 100 periods of the fundamental mode (nx = 200)")
print("=" * 78)
u0 = np.sin(np.pi * xg[1:-1]); v0 = np.zeros(m)
period = 2.0 / Cw
nu_ = 0.8
dt = nu_ * dx / Cw
nsteps = int(round(100 * period / dt))
print(f"dt = {dt:.5f}  (Courant {nu_}),  {nsteps:,} steps,  t_final = {nsteps*dt:.1f}\n")
print(f"{'scheme':<32} {'E_end/E_0':>12} {'drift':>12} {'osc. band':>12}")
E0 = energy(u0, v0)
res43 = {}
for name, fn in (("central difference (explicit)", lambda: run_newmark(u0, v0, dt, nsteps, 0.5, 0.0, 20)),
                 ("Newmark trapezoidal (1/2,1/4)", lambda: run_newmark(u0, v0, dt, nsteps, 0.5, 0.25, 20)),
                 ("Newmark damped (0.6,0.3025)", lambda: run_newmark(u0, v0, dt, nsteps, 0.6, 0.3025, 20)),
                 ("backward Euler", lambda: run_be(u0, v0, dt, nsteps, 20))):
    _, _, E = fn()
    res43[name] = E
    print(f"{name:<32} {E[-1]/E0:>12.6f} {E[-1]/E0 - 1:>12.2e} "
          f"{(E.max()-E.min())/E0:>12.2e}")
print("\nThe explicit central difference and trapezoidal Newmark hold energy to roundoff")
print("over ten thousand steps. They are symplectic: energy oscillates within a narrow")
print("band and never drifts. Backward Euler is unconditionally stable and bleeds the")
print("wave away entirely -- stability here means the answer decays to nothing.")

# =====================================================================
# 4.4  A pulse hitting a wall
# =====================================================================
print()
print("=" * 78)
print("4.4  Reflection: a Gaussian pulse striking a fixed and a free end")
print("=" * 78)
nx2 = 400
dx2 = 1.0 / nx2
x2 = np.linspace(0, 1, nx2 + 1)
pulse = np.exp(-((x2 - 0.5) / 0.04) ** 2)


def run_wave_bc(dt, nsteps, bc, scheme="leapfrog"):
    """Full-domain leapfrog with bc in {'fixed','free'} at BOTH ends."""
    u_old = pulse.copy()
    u = pulse.copy()                                  # v0 = 0: splits into two pulses
    nu2 = (Cw * dt / dx2) ** 2
    snaps = {}
    for n in range(nsteps + 1):
        if n in SNAP:
            snaps[n] = u.copy()
        lapu = np.zeros_like(u)
        lapu[1:-1] = u[2:] - 2 * u[1:-1] + u[:-2]
        if bc == "free":                              # ghost mirror: du/dx = 0
            lapu[0] = 2 * (u[1] - u[0]); lapu[-1] = 2 * (u[-2] - u[-1])
        u_new = 2 * u - u_old + nu2 * lapu
        if bc == "fixed":
            u_new[0] = 0.0; u_new[-1] = 0.0
        u_old, u = u, u_new
    return snaps


dt2 = 1.0 * dx2 / Cw                                  # Courant exactly 1
SNAP = {0, int(0.20 / dt2), int(0.45 / dt2), int(0.70 / dt2)}
snap_fixed = run_wave_bc(dt2, max(SNAP), "fixed")
snap_free = run_wave_bc(dt2, max(SNAP), "free")
print("Courant number exactly 1; the leapfrog scheme is then exact for this equation.")
for n in sorted(SNAP):
    t = n * dt2
    print(f"  t = {t:.3f}   fixed end: min u = {snap_fixed[n].min():+.4f}, "
          f"max u = {snap_fixed[n].max():+.4f}   |   free end: "
          f"min u = {snap_free[n].min():+.4f}, max u = {snap_free[n].max():+.4f}")
print("\nAt a fixed end the reflected pulse is INVERTED -- the wall must supply a force")
print("that cancels the incoming displacement. At a free end it is not. Both are")
print("captured exactly by an explicit scheme at Courant 1, and neither survives an")
print("implicit step large enough to be worth taking.")

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))

Om = np.linspace(0.01, 3.0, 400)
for (name, f), c in zip(SCH.items(), (BLU, GRN, PUR, ORA, RED)):
    pe = []
    for O in Om:
        rho, Omb = spectral(f(O))
        pe.append(100 * (O / Omb - 1) if rho <= 1 + 1e-9 else np.nan)
    ax[0].plot(Om / (2 * np.pi), pe, color=c, lw=1.3, label=name.split(" (")[0])
ax[0].axhline(0, color=GRY, lw=0.8)
ax[0].set_xlabel(r"$\Delta t / T$"); ax[0].set_ylabel("period error  (%)")
ax[0].set_title("Period elongation\n(curves stop where the scheme is unstable)", fontsize=9)
ax[0].legend(fontsize=7, frameon=False); ax[0].set_ylim(-15, 60); ax[0].grid(alpha=0.25)

tt = np.arange(len(res43["backward Euler"])) * 20 * dt / period
for (name, E), c in zip(res43.items(), (BLU, GRN, ORA, RED)):
    ax[1].plot(tt, E / E0, color=c, lw=1.2, label=name.split(" (")[0])
ax[1].set_xlabel("periods"); ax[1].set_ylabel("$E/E_0$")
ax[1].set_title("Energy over 100 periods", fontsize=9)
ax[1].legend(fontsize=7, frameon=False); ax[1].grid(alpha=0.25); ax[1].set_ylim(-0.05, 1.15)

for n, c, a_ in zip(sorted(SNAP), (GRY, BLU, ORA, RED), (0.5, 1, 1, 1)):
    ax[2].plot(x2, snap_fixed[n], color=c, lw=1.2, alpha=a_, label=f"$t$ = {n*dt2:.2f}")
ax[2].axhline(0, color=GRY, lw=0.7)
ax[2].set_xlabel("$x$"); ax[2].set_title("Pulse reflecting off fixed ends\n(inverts)", fontsize=9)
ax[2].legend(fontsize=7, frameon=False)

fig.tight_layout(); fig.savefig(f"{FIG}/time_fig5_wave.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigures -> {os.path.normpath(FIG)}")

# =====================================================================
# 4.1  The modes of the wave operator
# =====================================================================
print()
print("=" * 78)
print("4.1  Where the modes are")
print("=" * 78)
NXW = 16
dxw = 1.0 / NXW
Mw = NXW - 1
mw = np.arange(1, Mw + 1)
phiw = mw * np.pi * dxw
omega = (2 * Cw / dxw) * np.sin(phiw / 2)
print(f"nx = {NXW}: {Mw} spatial shapes, {2*Mw} eigenvalues (a +-i pair each)\n")
print(f"{'m':>4} {'phi':>8} {'omega_m':>10} {'lambda = +-i*omega':>22} "
      f"{'period 2pi/omega':>18}")
for m in (1, 2, 4, Mw - 1, Mw):
    print(f"{m:>4} {phiw[m-1]:>8.4f} {omega[m-1]:>10.2f} "
          f"{'+-' + str(round(omega[m-1],2)) + 'i':>22} "
          f"{2*np.pi/omega[m-1]:>18.4f}")
print(f"\n|lambda|_max = {omega[-1]:.2f} ~ 2c/dx = {2*Cw/dxw:.1f}   (ONE power of dx)")
print(f"|lambda|_min = {omega[0]:.2f} ~ c*pi = {Cw*np.pi:.2f}   (the fundamental)")
print(f"ratio = {omega[-1]/omega[0]:.1f}")
print("\nAll purely imaginary: nothing decays, ever. The exact |g| = 1 for every mode,")
print("so ANY amplitude error is a defect -- which is the opposite of section 2.")

fig, axw = plt.subplots(1, 2, figsize=(11.5, 3.9))
xw = np.linspace(0, 1, NXW + 1)
for m, c_, lbl in ((1, BLU, "$m=1$  fundamental"), (2, GRN, "$m=2$"),
                   (Mw, RED, f"$m={Mw}$  fastest")):
    v = np.zeros(NXW + 1); v[1:-1] = np.sin(m * np.pi * xw[1:-1])
    axw[0].plot(xw, v / np.max(np.abs(v)), "-o", ms=4, color=c_, lw=1.4,
                label=f"{lbl},  $\\omega$ = {omega[m-1]:.1f}")
axw[0].axhline(0, color=GRY, lw=0.7); axw[0].set_ylim(-1.15, 1.8)
axw[0].set_xlabel("$x$"); axw[0].set_ylabel("mode shape")
axw[0].set_title(f"Three of the {Mw} spatial shapes, $n_x$ = {NXW}", fontsize=10)
axw[0].legend(fontsize=7.5, frameon=False, loc="upper center")

axw[1].plot(np.zeros(Mw), omega, "o", ms=6, color=PUR)
axw[1].plot(np.zeros(Mw), -omega, "o", ms=6, color=PUR)
axw[1].axhline(0, color=GRY, lw=0.6); axw[1].axvline(0, color=GRY, lw=0.6)
axw[1].set_xlim(-12, 12); axw[1].set_ylim(-40, 40)
axw[1].set_xlabel(r"Re$(\lambda)$"); axw[1].set_ylabel(r"Im$(\lambda)$")
axw[1].set_title(r"All $2M$ eigenvalues: on the imaginary axis, in $\pm$ pairs",
                 fontsize=9.5)
axw[1].annotate("every eigenvalue has\nRe$(\\lambda) = 0$ exactly:\nno decay anywhere",
                xy=(0, 20), xytext=(-11, 34), fontsize=7.5, color=PUR, va="top",
                arrowprops=dict(arrowstyle="->", color=PUR, lw=0.7))
fig.suptitle("The wave equation: two eigenvalues per shape, none of them damped",
             y=1.04, fontsize=10)
fig.tight_layout()
fig.savefig(f"{FIG}/time_fig12_wave_modes.png", bbox_inches="tight")
plt.close(fig)

# =====================================================================
# 4.5  What accuracy needs — and it tightens with run length
# =====================================================================
print()
print("=" * 78)
print("4.5  Accuracy for waves is a PHASE budget, and it depends on how long you run")
print("=" * 78)
print("A period error of eps per cycle accumulates: after N cycles the wave is")
print("eps*N of a period out of phase. Requiring that to stay under 10%:\n")
print(f"{'scheme':<32} {'dt/T for 10% phase error after':>34}")
print(f"{'':32} {'1 cycle':>11} {'10':>11} {'100':>11}")
for name, f in SCH.items():
    row = f"{name:<32} "
    for N in (1, 10, 100):
        tol = 0.10 / N
        lo, hi = 1e-5, 0.49
        for _ in range(80):
            mid = np.sqrt(lo * hi)
            Om = 2 * np.pi * mid
            rho, Omb = spectral(f(Om))
            ok = rho <= 1 + 1e-9 and abs(Om / Omb - 1) <= tol
            lo, hi = (mid, hi) if ok else (lo, mid)
        row += f"{lo:>11.4f} "
    print(row)
print("\nRunning 100 times longer forces a step roughly 10x smaller for the second-order")
print("schemes (error ~ dt^2, so dt ~ sqrt(tol)). NOTHING like this happens in sec.2,")
print("where the solution decays and old errors decay with it.")

# =====================================================================
# 4.6  Recommended step across grids
# =====================================================================
print()
print("=" * 78)
print("4.6  Recommended step across grids (fundamental mode, 10 periods, rho = 4.5)")
print("=" * 78)
print(f"{'nx':>5} | {'scheme':<30} {'dt':>10} {'steps':>8} {'work':>9} {'E_end/E_0':>11}")
print("-" * 80)
for nxw2 in (50, 100, 200, 400):
    dxq = 1.0 / nxw2
    mq = nxw2 - 1
    xq = np.linspace(0, 1, nxw2 + 1)
    om_max = (2 * Cw / dxq) * np.sin((mq * np.pi * dxq) / 2)
    dt_stab = 2.0 / om_max
    Tend = 10 * 2.0 / Cw
    for label, gamma, beta, implicit in (
            ("central difference (explicit)", 0.5, 0.0, False),
            ("Newmark trapezoidal", 0.5, 0.25, True),
            ("Newmark trapezoidal, 4x step", 0.5, 0.25, True),
            ("backward Euler", None, None, True)):
        dt = dt_stab * 0.9 * (4 if "4x" in label else 1)
        n = max(1, int(round(Tend / dt)))
        globals()['nx'] = nxw2; globals()['dx'] = dxq; globals()['m'] = mq
        globals()['xg'] = xq
        u0q = np.sin(np.pi * xq[1:-1]); v0q = np.zeros(mq)
        E0q = energy(u0q, v0q)
        if gamma is None:
            _, _, E = run_be(u0q, v0q, Tend / n, n, max(1, n // 50))
        else:
            _, _, E = run_newmark(u0q, v0q, Tend / n, n, gamma, beta, max(1, n // 50))
        w = n * (RHO_W if implicit else 1.0)
        print(f"{nxw2 if 'central' in label else '':>5} | {label:<30} {Tend/n:>10.2e} "
              f"{n:>8} {w:>9.0f} {E[-1]/E0q:>11.6f}")
    print("-" * 80)
