import os, json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

# figures are written next to the document, as ../figs relative to this script
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figs")
_data = lambda name: os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "figure.dpi": 160,
                     "axes.spines.top": False, "axes.spines.right": False})

BLU, ORA, GRN, RED, GRY = "#2b6cb0", "#dd6b20", "#2f855a", "#c53030", "#718096"

# ----------------------------------------------------------------- fig 1
n = 4
h = 1.0 / n
fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.7))

for a in ax:
    a.set_xlim(-0.12, 1.12); a.set_ylim(-0.12, 1.12)
    a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values():
        s.set_visible(False)

# -- FD: nodes
a = ax[0]
for i in range(n + 1):
    a.plot([i * h, i * h], [0, 1], color=GRY, lw=0.6, zorder=1)
    a.plot([0, 1], [i * h, i * h], color=GRY, lw=0.6, zorder=1)
for i in range(n + 1):
    for j in range(n + 1):
        interior = 0 < i < n and 0 < j < n
        a.plot(i * h, j * h, "o", ms=13 if interior else 5,
               color=BLU if interior else "white",
               mec=BLU, mew=1.2, zorder=3)
        if interior:                       # dof index, ordered as in disc.py
            a.text(i * h, j * h, str((j - 1) * (n - 1) + (i - 1)), color="white",
                   fontsize=6.5, ha="center", va="center", zorder=5)
# highlight the stencil
ci = cj = 2
for di, dj in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)):
    a.plot((ci + di) * h, (cj + dj) * h, "o", ms=18, mfc="none", mec=RED, mew=1.6, zorder=4)
a.plot([(ci - 1) * h, (ci + 1) * h], [cj * h, cj * h], color=RED, lw=1.4, zorder=2)
a.plot([ci * h, ci * h], [(cj - 1) * h, (cj + 1) * h], color=RED, lw=1.4, zorder=2)
a.set_title("Finite difference\nunknowns at nodes  •  9 interior dof", color=BLU)

# -- FVM: cells
a = ax[1]
for i in range(n):
    for j in range(n):
        a.add_patch(Rectangle((i * h, j * h), h, h, fc="none", ec=GRY, lw=0.6, zorder=1))
        a.plot((i + .5) * h, (j + .5) * h, "s", ms=14, color=ORA, zorder=3)
        a.text((i + .5) * h, (j + .5) * h, str(j * n + i), color="white",
               fontsize=6.5, ha="center", va="center", zorder=5)
ci = cj = 1
a.add_patch(Rectangle((ci * h, cj * h), h, h, fc=ORA, alpha=0.16, ec=ORA, lw=1.8, zorder=2))
for (dx, dy, ang) in ((h, .5 * h, 0), (0, .5 * h, 180), (.5 * h, h, 90), (.5 * h, 0, 270)):
    x0, y0 = ci * h + dx, cj * h + dy
    dxx, dyy = 0.055 * np.cos(np.deg2rad(ang)), 0.055 * np.sin(np.deg2rad(ang))
    a.arrow(x0, y0, dxx, dyy, head_width=0.028, color=RED, lw=1.3, zorder=5,
            length_includes_head=True)
a.set_title("Finite volume\nunknowns are cell averages  •  16 dof", color=ORA)

# -- FEM: P1 triangles
a = ax[2]
for i in range(n):
    for j in range(n):
        p = [(i * h, j * h), ((i + 1) * h, j * h), ((i + 1) * h, (j + 1) * h), (i * h, (j + 1) * h)]
        a.add_patch(Polygon([p[0], p[1], p[2]], fc="none", ec=GRY, lw=0.6, zorder=1))
        a.add_patch(Polygon([p[0], p[2], p[3]], fc="none", ec=GRY, lw=0.6, zorder=1))
for i in range(n + 1):
    for j in range(n + 1):
        interior = 0 < i < n and 0 < j < n
        a.plot(i * h, j * h, "^", ms=15 if interior else 5,
               color=GRN if interior else "white", mec=GRN, mew=1.2, zorder=3)
        if interior:
            a.text(i * h, j * h - 0.06 * h, str((j - 1) * (n - 1) + (i - 1)), color="white",
                   fontsize=6.5, ha="center", va="center", zorder=5)
ci = cj = 2
# the six triangles meeting at the centre node: E, NE, N, W, SW, S.
# the centre itself is interior to the patch and is not a boundary vertex.
patch = [((ci + 1) * h, cj * h), ((ci + 1) * h, (cj + 1) * h), (ci * h, (cj + 1) * h),
         ((ci - 1) * h, cj * h), ((ci - 1) * h, (cj - 1) * h), (ci * h, (cj - 1) * h)]
a.add_patch(Polygon(patch, fc=GRN, alpha=0.15, zorder=2))
a.add_patch(Polygon(patch, fc="none", ec=RED, lw=1.6, zorder=4))
a.set_title("Finite element (P1)\nunknowns at nodes  •  9 interior dof", color=GRN)

fig.suptitle("The same 4 x 4 discretization of the unit square, three ways", y=1.0)
fig.text(0.5, -0.03, "red: the support of one equation — a stencil, a control volume, "
                     "an element patch.\nNumbers are the unknown ordering used in section 2, "
                     "i running fastest", ha="center", color=RED, fontsize=8.5)
fig.tight_layout()
fig.savefig(f"{OUT}/fig1_three_meshes.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------- fig 2
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.5))
s5 = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]], float)
s9 = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]], float) / 3.0
for a, S, ttl, sub in ((ax[0], s5, "FD  =  FVM  =  FEM P1",
                        "5-point,  identical to machine precision"),
                       (ax[1], s9, "FEM Q1 (bilinear quads)",
                        "9-point,  entries $\\times \\frac{1}{3}$")):
    a.imshow(np.abs(S), cmap="Blues", vmin=0, vmax=3)
    for i in range(3):
        for j in range(3):
            v = S[i, j]
            txt = "0" if abs(v) < 1e-12 else (f"{v:.0f}" if abs(v - round(v)) < 1e-9
                                              else f"{v:+.3f}")
            a.text(j, i, txt, ha="center", va="center", fontsize=13,
                   color="white" if abs(v) > 1.5 else "#1a202c",
                   fontweight="bold" if i == j == 1 else "normal")
    a.set_xticks([]); a.set_yticks([])
    a.set_title(ttl + "\n" + sub, fontsize=9)
fig.suptitle("Interior stencil for $-\\nabla^2$ on a uniform mesh, scaled by $h^2$", y=1.02)
fig.tight_layout()
fig.savefig(f"{OUT}/fig2_stencils.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------- fig 3
conv = json.load(open(_data("conv.json")))
fig, ax = plt.subplots(figsize=(5.4, 4.2))
style = {"FD  5-point": (BLU, "o", "-"), "FVM cell-centred": (ORA, "s", "-"),
         "FEM P1 triangles": (GRN, "^", "-"), "FEM Q1 quads": ("#805ad5", "d", "--")}
for k, v in conv.items():
    c, m, ls = style[k]
    ax.loglog(v["h"], v["L2"], ls, marker=m, color=c, label=k, ms=5, lw=1.4)
hh = np.array(conv["FD  5-point"]["h"])
ax.loglog(hh, 0.42 * hh**2, ":", color=GRY, lw=1.2)
ax.text(hh[1], 0.42 * hh[1]**2 * 1.5, r"slope 2", color=GRY, fontsize=9)
ax.set_xlabel("$h$"); ax.set_ylabel(r"$L^2$ error")
ax.set_title("All four are second order on the smooth problem")
ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.25, which="both")
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_convergence.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------- fig 4
adv = json.load(open(_data("adv.json")))
Pe = 50.0
xf = np.linspace(0, 1, 400)
Tex = (np.exp(Pe * xf) - 1) / (np.exp(Pe) - 1)
fig, ax = plt.subplots(1, 2, figsize=(9.6, 3.9), sharey=True)
for a in ax:
    a.plot(xf, Tex, color="k", lw=1.6, label="exact")
    a.axhline(0, color=GRY, lw=0.6)
    a.set_xlabel("$x$"); a.grid(alpha=0.22)
ax[0].plot(adv["x_c"], adv["central"], "o--", color=RED, ms=5, label="FVM central")
ax[0].plot(adv["x_c"], adv["upwind"], "s-", color=ORA, ms=5, label="FVM upwind")
ax[0].set_title("Finite volume")
ax[1].plot(adv["x_n"], adv["galerkin"], "o--", color=RED, ms=5, label="FEM Galerkin")
ax[1].plot(adv["x_n"], adv["supg"], "^-", color=GRN, ms=5, label="FEM SUPG")
ax[1].set_title("Finite element")
ax[0].set_ylabel("$T$")
for a in ax:
    a.legend(fontsize=8, frameon=False, loc="upper left")
fig.suptitle(f"Advection–diffusion, global Pe = 50, 10 cells, cell Pe$_h$ = "
             f"{adv['Pe_h']:.0f}", y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_advection.png", bbox_inches="tight")
plt.close(fig)

# ----------------------------------------------------------------- fig 5
sl = json.load(open(_data("slab.json")))
q = sl["q_exact"]; K1, K2 = sl["K1"], sl["K2"]
xe = np.linspace(0, 1, 400)
Te = np.where(xe < 0.5, q * xe / K1, q * 0.5 / K1 + q * (xe - 0.5) / K2)
fig, ax = plt.subplots(figsize=(5.8, 4.0))
ax.axvspan(0, 0.5, color=BLU, alpha=0.06)
ax.axvspan(0.5, 1, color=ORA, alpha=0.08)
ax.plot(xe, Te, "k", lw=1.6, label="exact")
ax.plot(sl["xc"], sl["harm"], "o-", color=GRN, ms=6, label="FVM, harmonic face $k$")
ax.plot(sl["xc"], sl["arith"], "s--", color=RED, ms=6, label="FVM, arithmetic face $k$")
ax.axvline(0.5, color=GRY, lw=1.0, ls=":")
ax.text(0.22, 0.15, f"$k = {K1:g}$", color=BLU, fontsize=10)
ax.text(0.72, 0.15, f"$k = {K2:g}$", color=ORA, fontsize=10)
ax.set_xlabel("$x$"); ax.set_ylabel("$T$")
ax.set_title("Two-material slab, 8 cells\nflux error: harmonic $9\\times10^{-16}$, "
             "arithmetic $14.2\\%$", fontsize=9.5)
ax.legend(fontsize=8, frameon=False, loc="lower right"); ax.grid(alpha=0.22)
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_interface.png", bbox_inches="tight")
plt.close(fig)

print("figures written to", OUT)
print("\n".join(sorted(os.listdir(OUT))))
