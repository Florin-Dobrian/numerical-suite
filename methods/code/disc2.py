"""E3 (discontinuous conductivity) and E4 (advection-diffusion) for DISCRETIZATION.md."""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_data = lambda name: os.path.join(_HERE, name)
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import json

np.set_printoptions(precision=6, suppress=True, linewidth=140)

# ======================================================================
# E3.  Two-material slab.  k = k1 for x<0.5, k2 for x>0.5.
#      T=0 at x=0, T=1 at x=1, insulated at y=0 and y=1.
#      Exact: piecewise linear, uniform flux q = 1/(0.5/k1 + 0.5/k2).
# ======================================================================
K1, K2 = 1.0, 1000.0
q_exact = 1.0 / (0.5 / K1 + 0.5 / K2)
T_iface_exact = q_exact * 0.5 / K1

print("=" * 78)
print("E3.  Discontinuous conductivity, k1 = %g, k2 = %g" % (K1, K2))
print("=" * 78)
print("exact flux through any vertical plane   q = %.9f" % q_exact)
print("exact temperature at the interface         = %.9f" % T_iface_exact)


def fv_slab(n, face_avg):
    """Cell-centred FVM, n x n cells, interface lands exactly on a cell face."""
    h = 1.0 / n
    xc = (np.arange(n) + 0.5) * h
    kc = np.where(xc < 0.5, K1, K2)
    idx = lambda i, j: j * n + i
    A = sp.lil_matrix((n * n, n * n))
    b = np.zeros(n * n)
    kf = (lambda a, c: 2 * a * c / (a + c)) if face_avg == "harmonic" \
        else (lambda a, c: 0.5 * (a + c))
    for j in range(n):
        for i in range(n):
            p = idx(i, j)
            diag = 0.0
            # x-direction
            if i > 0:
                t = kf(kc[i], kc[i - 1]); A[p, idx(i - 1, j)] -= t; diag += t
            else:
                t = 2 * kc[i]; diag += t                       # T=0 wall
            if i < n - 1:
                t = kf(kc[i], kc[i + 1]); A[p, idx(i + 1, j)] -= t; diag += t
            else:
                t = 2 * kc[i]; diag += t; b[p] += t * 1.0      # T=1 wall
            # y-direction: interior faces only (insulated top and bottom)
            for jj in (j - 1, j + 1):
                if 0 <= jj < n:
                    t = kf(kc[i], kc[i]); A[p, idx(i, jj)] -= t; diag += t
            A[p, p] = diag
    T = spla.spsolve(A.tocsr(), b).reshape(n, n)
    # flux through the wall at x=0, per unit height:  sum over j of k*(T_P-0)/(h/2)*h
    q_wall = np.sum(2 * kc[0] * (T[:, 0] - 0.0))
    # flux through the material interface face (between cell n/2-1 and n/2)
    a = n // 2 - 1
    q_if = np.sum(kf(kc[a], kc[a + 1]) * (T[:, a] - T[:, a + 1]))
    return T, xc, q_wall, q_if


for n in (8, 16, 32):
    Th, xc, qw_h, qi_h = fv_slab(n, "harmonic")
    Ta, _, qw_a, qi_a = fv_slab(n, "arithmetic")
    print(f"\nn = {n:3d}   (h = {1/n:.4f})")
    print(f"  harmonic  face k :  q at wall = {qw_h:.9f}   "
          f"rel err = {abs(qw_h-q_exact)/q_exact:.3e}")
    print(f"  arithmetic face k:  q at wall = {qw_a:.9f}   "
          f"rel err = {abs(qw_a-q_exact)/q_exact:.3e}")
    print(f"  harmonic  : T just left of interface = {Th[0, n//2-1]:.6f}")
    print(f"  arithmetic: T just left of interface = {Ta[0, n//2-1]:.6f}")

# local conservation check on the harmonic FVM solution: is the flux the same
# through every vertical plane?
Th, xc, _, _ = fv_slab(32, "harmonic")
n = 32
fluxes = []
kc = np.where(xc < 0.5, K1, K2)
kfh = lambda a, c: 2 * a * c / (a + c)
for a in range(n - 1):
    fluxes.append(np.sum(kfh(kc[a], kc[a + 1]) * (Th[:, a] - Th[:, a + 1])))
fluxes = np.array(fluxes)
print("\nFVM, harmonic, n=32: flux through each of the 31 interior vertical planes")
print("  min = %.12f   max = %.12f   spread = %.3e"
      % (fluxes.min(), fluxes.max(), fluxes.max() - fluxes.min()))
print("  (a uniform flux is the exact answer; the spread is the conservation defect)")


def p1_slab(n):
    """FEM P1 on an interface-aligned triangulation, element-wise constant k."""
    h = 1.0 / n
    nid = lambda i, j: j * (n + 1) + i
    nn = (n + 1)**2
    A = sp.lil_matrix((nn, nn))
    F = np.zeros(nn)
    for j in range(n):
        for i in range(n):
            kel = K1 if (i + 0.5) * h < 0.5 else K2
            c = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            for tri in ([c[0], c[1], c[2]], [c[0], c[2], c[3]]):
                g = [nid(*t) for t in tri]
                xs = np.array([t[0] * h for t in tri])
                ys = np.array([t[1] * h for t in tri])
                bb = np.array([ys[1] - ys[2], ys[2] - ys[0], ys[0] - ys[1]])
                cc = np.array([xs[2] - xs[1], xs[0] - xs[2], xs[1] - xs[0]])
                area = 0.5 * abs(bb[0] * cc[1] - bb[1] * cc[0])
                Ke = kel * (np.outer(bb, bb) + np.outer(cc, cc)) / (4 * area)
                for a in range(3):
                    for d in range(3):
                        A[g[a], g[d]] += Ke[a, d]
    A = A.tocsr()
    dirich = {}
    for j in range(n + 1):
        dirich[nid(0, j)] = 0.0
        dirich[nid(n, j)] = 1.0
    Afull = A.copy()
    for p, v in dirich.items():
        F -= np.asarray(Afull[:, p].todense()).ravel() * v
    A = A.tolil()
    for p, v in dirich.items():
        A.rows[p] = [p]; A.data[p] = [1.0]; F[p] = v
    A = A.tocsr().tolil()
    for p in dirich:                      # zero the columns too, then restore diag
        for r in range(A.shape[0]):
            if r != p and p in A.rows[r]:
                k = A.rows[r].index(p); A.rows[r].pop(k); A.data[r].pop(k)
    T = spla.spsolve(A.tocsr(), F)
    # reaction flux at the x=1 wall = sum of residual entries there
    R = Afull @ T
    q = -np.sum(R[[nid(n, j) for j in range(n + 1)]])
    return T.reshape(n + 1, n + 1), q


for n in (8, 16, 32):
    T, q = p1_slab(n)
    print(f"FEM P1, interface-aligned, n = {n:3d}: reaction flux at x=1 = {abs(q):.9f}"
          f"   rel err = {abs(abs(q)-q_exact)/q_exact:.3e}")

# ======================================================================
# E4.  Advection-diffusion.  u dT/dx = alpha d2T/dx2, T(0)=0, T(1)=1.
#      2D domain, insulated top and bottom, so the solution is x-only.
# ======================================================================
print()
print("=" * 78)
print("E4.  Advection-diffusion at high cell Peclet number")
print("=" * 78)
U, L = 1.0, 1.0
Pe_global = 50.0
alpha = U * L / Pe_global
Texact = lambda x: (np.exp(Pe_global * x) - 1) / (np.exp(Pe_global) - 1)


def fv_ad(n, scheme):
    """Cell-centred FVM, assembled face by face from the outward flux balance.
    Walls are Dirichlet: the face there takes the wall value, and the diffusive
    distance is h/2. The 2-D problem is y-invariant, so this is its x-line."""
    h = 1.0 / n
    D = alpha / h
    A = np.zeros((n, n)); b = np.zeros(n)
    for i in range(n):
        # ---- east face ----
        if i < n - 1:
            A[i, i] += D; A[i, i + 1] -= D                    # -De(T_E - T_P)
            if scheme == "upwind":
                A[i, i] += U                                  # +U T_P
            else:
                A[i, i] += 0.5 * U; A[i, i + 1] += 0.5 * U
        else:
            A[i, i] += 2 * D; b[i] += 2 * D * 1.0             # wall T = 1
            A[i, i] += U                     # outflow: the face carries T_P
        # ---- west face ----
        if i > 0:
            A[i, i] += D; A[i, i - 1] -= D                    # -Dw(T_W - T_P)
            if scheme == "upwind":
                A[i, i - 1] -= U                              # -U T_W
            else:
                A[i, i] -= 0.5 * U; A[i, i - 1] -= 0.5 * U
        else:
            A[i, i] += 2 * D; b[i] += 2 * D * 0.0             # wall T = 0
            b[i] += U * 0.0                  # inflow: the face carries T_wall = 0
    return np.linalg.solve(A, b), (np.arange(n) + 0.5) * h


def p1_ad(n, supg=False):
    """P1 Galerkin (= central) and SUPG, 1-D-in-x."""
    h = 1.0 / n
    nn = n + 1
    A = sp.lil_matrix((nn, nn)); F = np.zeros(nn)
    Pe_h = U * h / (2 * alpha)
    tau = (h / (2 * U)) * (1 / np.tanh(Pe_h) - 1 / Pe_h) if supg else 0.0
    for e in range(n):
        g = [e, e + 1]
        Kd = alpha / h * np.array([[1, -1], [-1, 1]])
        Ka = U * np.array([[-0.5, 0.5], [-0.5, 0.5]])
        Ks = tau * U * U / h * np.array([[1, -1], [-1, 1]]) if supg else np.zeros((2, 2))
        Ke = Kd + Ka + Ks
        for a in range(2):
            for d in range(2):
                A[g[a], g[d]] += Ke[a, d]
    Afull = A.tocsr().copy()
    A = A.tolil()
    for p, v in ((0, 0.0), (nn - 1, 1.0)):
        F -= np.asarray(Afull[:, p].todense()).ravel() * v
    for p, v in ((0, 0.0), (nn - 1, 1.0)):
        A.rows[p] = [p]; A.data[p] = [1.0]; F[p] = v
    return spla.spsolve(A.tocsr(), F), np.arange(nn) * h


out = {}
for n in (10, 20, 50):
    Pe_h = U * (1.0 / n) / alpha
    Tc, xc = fv_ad(n, "central")
    Tu, _ = fv_ad(n, "upwind")
    Tg, xn = p1_ad(n, supg=False)
    Ts, _ = p1_ad(n, supg=True)
    print(f"\nn = {n:3d}  cell Peclet = U h / alpha = {Pe_h:.2f}")
    print(f"  FVM central   : min = {Tc.min():9.4f}  max = {Tc.max():9.4f}  "
          f"{'OSCILLATES' if Tc.min() < -1e-8 else 'monotone'}")
    print(f"  FVM upwind    : min = {Tu.min():9.4f}  max = {Tu.max():9.4f}  "
          f"{'OSCILLATES' if Tu.min() < -1e-8 else 'monotone'}")
    print(f"  FEM Galerkin  : min = {Tg.min():9.4f}  max = {Tg.max():9.4f}  "
          f"{'OSCILLATES' if Tg.min() < -1e-8 else 'monotone'}")
    print(f"  FEM SUPG      : min = {Ts.min():9.4f}  max = {Ts.max():9.4f}  "
          f"{'OSCILLATES' if Ts.min() < -1e-8 else 'monotone'}")
    print(f"  upwind error at the last interior node: "
          f"{abs(Tu[-1]-Texact(xc[-1])):.4e}   (numerical diffusion)")
    if n == 10:
        out = {'x_c': xc.tolist(), 'x_n': xn.tolist(),
               'central': Tc.tolist(), 'upwind': Tu.tolist(),
               'galerkin': Tg.tolist(), 'supg': Ts.tolist(),
               'Pe_h': Pe_h}
json.dump(out, open(_data('adv.json'), 'w'))
print("\nSUPG tau uses the exact-nodal (doubly asymptotic) formula "
      "tau = (h/2U)(coth(Pe_h/2) - 2/Pe_h).")

print("\nWhere exactly does the central scheme lose monotonicity?")
for n in (24, 25, 26, 30):
    Pe_h = U * (1.0 / n) / alpha
    Tc, _ = fv_ad(n, "central")
    print(f"  n = {n:3d}  Pe_h = {Pe_h:5.3f}   min(T) = {Tc.min():10.3e}  "
          f"{'oscillates' if Tc.min() < -1e-12 else 'monotone'}")

Th8, xc8, _, _ = fv_slab(8, "harmonic")
Ta8, _, _, _ = fv_slab(8, "arithmetic")
json.dump({'xc': xc8.tolist(), 'harm': Th8[0].tolist(), 'arith': Ta8[0].tolist(),
           'q_exact': q_exact, 'K1': K1, 'K2': K2}, open(_data('slab.json'), 'w'))
