"""
TIME_INTEGRATION.md sec.1.4 -- stability and accuracy are different questions.

One step multiplies a Fourier mode by g. Stability asks only  |g| <= 1.
Accuracy asks whether g is close to the exact factor e^z.  A scheme can pass
either test while failing the other, and the pure oscillation lambda = i*omega
separates them completely.
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

R = {"Forward Euler":  lambda z: 1 + z,
     "Backward Euler": lambda z: 1 / (1 - z),
     "Crank-Nicolson": lambda z: (1 + z / 2) / (1 - z / 2)}

# ---------------------------------------------------------------- 1. decay
print("=" * 78)
print("A.  A DECAYING mode, z = lambda dt = -2.  Exact factor e^-2 = 0.135335")
print("=" * 78)
print(f"{'scheme':<17} {'g':>12} {'|g|':>8} {'stable?':>9} {'|g - e^z|':>12}   accuracy")
ez = np.exp(-2)
for n, f in R.items():
    g = f(-2 + 0j).real
    print(f"{n:<17} {g:>12.6f} {abs(g):>8.4f} {str(abs(g) <= 1):>9} {abs(g-ez):>12.6f}"
          f"   {'off by %.1fx' % (g/ez) if g > 0 else 'WRONG SIGN'}")
print("\nAll three are stable at z = -2. None is accurate. Backward Euler overestimates")
print("the surviving amplitude by 2.5x, Crank-Nicolson annihilates it, forward Euler")
print("returns the negative of it. Stability said nothing about any of that.")

# ---------------------------------------------------------------- 2. oscillation
print()
print("=" * 78)
print("B.  A PURE OSCILLATION, z = i*Omega.  The exact factor has |e^z| = 1 EXACTLY")
print("=" * 78)
print("Nothing decays: the true solution keeps its amplitude forever. So |g| < 1 is")
print("as wrong as |g| > 1 -- it is just wrong in the direction stability rewards.\n")
print(f"{'Omega':>7}  " + "".join(f"{n:>34}" for n in R))
print(f"{'':>7}  " + "".join(f"{'|g|':>12}{'amp/cycle':>12}{'verdict':>10}" for _ in R))
for Om in (0.1, 0.5, 1.0, 2.0):
    row = f"{Om:>7.2f}  "
    for n, f in R.items():
        g = f(1j * Om)
        ncyc = 2 * np.pi / Om
        amp = abs(g) ** ncyc
        v = "grows" if abs(g) > 1 + 1e-12 else ("decays" if abs(g) < 1 - 1e-12 else "exact")
        row += f"{abs(g):>12.6f}{amp:>12.4f}{v:>10}"
    print(row)
print("\nForward Euler  |g| = sqrt(1 + Omega^2) > 1 at EVERY step size: unstable always.")
print("Backward Euler |g| = 1/sqrt(1 + Omega^2) < 1 at every step size: stable always,")
print("               and it eats the oscillation it was asked to compute.")
print("Crank-Nicolson |g| = 1 identically: neutrally stable, amplitude exactly right.")
print("               Its error is entirely in the PHASE -- and that is second order.")

print()
print("=" * 78)
print("C.  Where accuracy goes when stability is guaranteed")
print("=" * 78)
print("Backward Euler on z = i*Omega. It never diverges. Watch the answer anyway:\n")
print(f"{'Omega':>7} {'steps/cycle':>12} {'|g|':>10} {'amplitude left after 10 cycles':>32}")
for Om in (0.05, 0.1, 0.25, 0.5, 1.0, 2.0):
    g = abs(R["Backward Euler"](1j * Om))
    print(f"{Om:>7.2f} {2*np.pi/Om:>12.1f} {g:>10.6f} {g**(10*2*np.pi/Om):>32.6f}")
print("\nEvery row is stable. The last row has destroyed 99.99% of the signal.")
print("'Unconditionally stable' is a promise about divergence, not about the answer.")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
Om = np.linspace(0, 3, 400)

for (n, f), c in zip(R.items(), (RED, BLU, GRN)):
    ax[0].plot(Om, np.abs(f(1j * Om)), color=c, lw=1.5, label=n)
ax[0].axhline(1, color="k", ls=":", lw=1.2)
ax[0].text(2.05, 1.08, "exact: $|g| = 1$", fontsize=8)
ax[0].fill_between(Om, 1, 2.5, color=RED, alpha=0.06)
ax[0].text(1.75, 1.75, "unstable region", color=RED, fontsize=8.5)
ax[0].text(0.35, 0.22, "stable — by damping a wave\nthat should never damp", color=BLU, fontsize=8)
ax[0].set_xlabel(r"$\Omega = \omega\Delta t$"); ax[0].set_ylabel("$|g|$")
ax[0].set_title("STABILITY: is $|g|\\leq1$?", fontsize=9.5)
ax[0].set_ylim(0, 2.5); ax[0].legend(fontsize=8, frameon=False, loc="lower left",
                                    bbox_to_anchor=(0.0, 0.62))

for (n, f), c in zip(R.items(), (RED, BLU, GRN)):
    ax[1].loglog(Om[1:], np.abs(f(1j * Om[1:]) - np.exp(1j * Om[1:])), color=c, lw=1.5, label=n)
ax[1].loglog(Om[1:], 0.5 * Om[1:]**2, ":", color=GRY, lw=1.1)
ax[1].loglog(Om[1:], 0.08 * Om[1:]**3, "--", color=GRY, lw=1.1)
ax[1].text(0.5, 0.09, "slope 2", color=GRY, fontsize=8)
ax[1].text(1.1, 0.02, "slope 3", color=GRY, fontsize=8)
ax[1].set_xlabel(r"$\Omega = \omega\Delta t$"); ax[1].set_ylabel(r"$|g - e^{i\Omega}|$")
ax[1].set_title("ACCURACY: is $g$ close to $e^{z}$?", fontsize=9.5)
ax[1].legend(fontsize=8, frameon=False); ax[1].grid(alpha=0.25, which="both")

ncyc = 10 * 2 * np.pi / Om[1:]
for (n, f), c in zip(R.items(), (RED, BLU, GRN)):
    ax[2].semilogy(Om[1:], np.clip(np.abs(f(1j * Om[1:]))**ncyc, 1e-8, 1e8),
                   color=c, lw=1.5, label=n)
ax[2].axhline(1, color="k", ls=":", lw=1.2)
ax[2].set_xlabel(r"$\Omega = \omega\Delta t$")
ax[2].set_ylabel("amplitude after 10 cycles")
ax[2].set_title("What that costs over 10 cycles", fontsize=9.5)
ax[2].set_ylim(1e-8, 1e8); ax[2].legend(fontsize=8, frameon=False); ax[2].grid(alpha=0.25)

fig.suptitle("Pure oscillation $\\lambda = i\\omega$: the two questions come apart completely",
             y=1.04, fontsize=10)
fig.tight_layout(); fig.savefig(f"{FIG}/time_fig7_stability_vs_accuracy.png", bbox_inches="tight")
plt.close(fig)
print(f"\nfigure -> {os.path.normpath(FIG)}/time_fig7_stability_vs_accuracy.png")
