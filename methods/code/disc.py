"""
Experiments for DISCRETIZATION.md.

One 2D domain, one PDE, three discretizations. Everything below is computed,
not asserted; numbers quoted in the document come from this script.
"""
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
_data = lambda name: os.path.join(_HERE, name)
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

np.set_printoptions(precision=6, suppress=True, linewidth=140)

# ----------------------------------------------------------------------
# Manufactured solution on the unit square:  -Lap u = f,  u = 0 on dOmega
#   u = sin(pi x) sin(pi y),  f = 2 pi^2 sin(pi x) sin(pi y)
# ----------------------------------------------------------------------
uex = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
fsrc = lambda x, y: 2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)

# exact cell average of u over [a,b]x[c,d], for the FVM comparison
def uex_cellavg(xl, xr, yl, yr):
    ix = (np.cos(np.pi * xl) - np.cos(np.pi * xr)) / np.pi
    iy = (np.cos(np.pi * yl) - np.cos(np.pi * yr)) / np.pi
    return ix * iy / ((xr - xl) * (yr - yl))


# ======================================================================
# 1. FINITE DIFFERENCE  -- 5-point, node-based, interior nodes unknown
# ======================================================================
def fd_solve(n):
    h = 1.0 / n
    m = n - 1                                   # interior nodes per direction
    idx = lambda i, j: (j - 1) * m + (i - 1)
    A = sp.lil_matrix((m * m, m * m))
    b = np.zeros(m * m)
    for j in range(1, n):
        for i in range(1, n):
            p = idx(i, j)
            A[p, p] = 4.0 / h**2
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 1 <= ii <= n - 1 and 1 <= jj <= n - 1:
                    A[p, idx(ii, jj)] = -1.0 / h**2
                # else: Dirichlet value 0, contributes nothing to b
            b[p] = fsrc(i * h, j * h)           # POINT value of f
    u = spla.spsolve(A.tocsr(), b)
    X, Y = np.meshgrid(np.arange(1, n) * h, np.arange(1, n) * h)
    return u, uex(X, Y).ravel(), h, A.tocsr()


# ======================================================================
# 2. FINITE VOLUME  -- cell-centred, n x n cells
# ======================================================================
def fv_solve(n, kfun=None, face_avg="harmonic"):
    h = 1.0 / n
    idx = lambda i, j: j * n + i
    A = sp.lil_matrix((n * n, n * n))
    b = np.zeros(n * n)
    xc = (np.arange(n) + 0.5) * h
    kc = np.ones((n, n)) if kfun is None else np.array(
        [[kfun(xc[i], xc[j]) for i in range(n)] for j in range(n)])

    def kface(k1, k2):
        if face_avg == "harmonic":
            return 2 * k1 * k2 / (k1 + k2)
        return 0.5 * (k1 + k2)

    for j in range(n):
        for i in range(n):
            p = idx(i, j)
            diag = 0.0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:
                    kf = kface(kc[j, i], kc[jj, ii])
                    t = kf                       # (k * face length h) / h
                    A[p, idx(ii, jj)] -= t
                    diag += t
                else:
                    # Dirichlet 0 at the wall: half-cell distance h/2
                    t = 2.0 * kc[j, i]
                    diag += t
            A[p, p] = diag
            # CELL AVERAGE of f, times cell volume h^2  ->  integral of f
            b[p] = h**2 * cellavg_f(i * h, (i + 1) * h, j * h, (j + 1) * h)
    u = spla.spsolve(A.tocsr(), b)
    X, Y = np.meshgrid(xc, xc)
    ubar = np.array([[uex_cellavg(i * h, (i + 1) * h, j * h, (j + 1) * h)
                      for i in range(n)] for j in range(n)]).ravel()
    return u, ubar, h, A.tocsr()


def cellavg_f(xl, xr, yl, yr):
    ix = (np.cos(np.pi * xl) - np.cos(np.pi * xr)) / np.pi
    iy = (np.cos(np.pi * yl) - np.cos(np.pi * yr)) / np.pi
    return 2 * np.pi**2 * ix * iy / ((xr - xl) * (yr - yl))


# ======================================================================
# 3a. FINITE ELEMENT -- P1 linear triangles, each square split by one diagonal
# ======================================================================
def p1_solve(n):
    h = 1.0 / n
    nn = (n + 1)**2
    nid = lambda i, j: j * (n + 1) + i
    A = sp.lil_matrix((nn, nn))
    F = np.zeros(nn)
    for j in range(n):
        for i in range(n):
            c = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
            for tri in ([c[0], c[1], c[2]], [c[0], c[2], c[3]]):
                g = [nid(*t) for t in tri]
                xs = np.array([t[0] * h for t in tri])
                ys = np.array([t[1] * h for t in tri])
                bcoef = np.array([ys[1] - ys[2], ys[2] - ys[0], ys[0] - ys[1]])
                ccoef = np.array([xs[2] - xs[1], xs[0] - xs[2], xs[1] - xs[0]])
                area = 0.5 * abs(bcoef[0] * ccoef[1] - bcoef[1] * ccoef[0])
                Ke = np.outer(bcoef, bcoef) + np.outer(ccoef, ccoef)
                Ke /= (4 * area)
                # load: 3-point midside rule (exact for quadratics)
                Fe = np.zeros(3)
                mids = [(0, 1), (1, 2), (2, 0)]
                for a, bb in mids:
                    xm, ym = 0.5 * (xs[a] + xs[bb]), 0.5 * (ys[a] + ys[bb])
                    w = area / 3.0
                    Fe[a] += w * fsrc(xm, ym) * 0.5
                    Fe[bb] += w * fsrc(xm, ym) * 0.5
                for a in range(3):
                    F[g[a]] += Fe[a]
                    for bb in range(3):
                        A[g[a], g[bb]] += Ke[a, bb]
    return _apply_dirichlet_and_solve(A, F, n, h)


# ======================================================================
# 3b. FINITE ELEMENT -- Q1 bilinear quadrilaterals, 2x2 Gauss
# ======================================================================
def q1_solve(n):
    h = 1.0 / n
    nn = (n + 1)**2
    nid = lambda i, j: j * (n + 1) + i
    A = sp.lil_matrix((nn, nn))
    F = np.zeros(nn)
    gp = np.array([-1, 1]) / np.sqrt(3.0)
    for j in range(n):
        for i in range(n):
            g = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)]
            x0, y0 = i * h, j * h
            Ke = np.zeros((4, 4))
            Fe = np.zeros(4)
            for xi in gp:
                for eta in gp:
                    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
                    dNxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
                    dNet = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
                    J = h / 2.0
                    B = np.vstack([dNxi / J, dNet / J])
                    detJ = J * J
                    Ke += (B.T @ B) * detJ
                    xg = x0 + h * (xi + 1) / 2
                    yg = y0 + h * (eta + 1) / 2
                    Fe += N * fsrc(xg, yg) * detJ
            for a in range(4):
                F[g[a]] += Fe[a]
                for bb in range(4):
                    A[g[a], g[bb]] += Ke[a, bb]
    return _apply_dirichlet_and_solve(A, F, n, h)


def _apply_dirichlet_and_solve(A, F, n, h):
    nn = (n + 1)**2
    nid = lambda i, j: j * (n + 1) + i
    bnd = set()
    for i in range(n + 1):
        bnd |= {nid(i, 0), nid(i, n), nid(0, i), nid(n, i)}
    Araw = A.tocsr().copy()
    A = A.tolil()
    for p in sorted(bnd):
        A.rows[p] = [p]
        A.data[p] = [1.0]
        F[p] = 0.0
    A = A.tocsr()
    for p in sorted(bnd):                      # symmetrise (zero BC: no lifting)
        pass
    u = spla.spsolve(A, F)
    free = np.array(sorted(set(range(nn)) - bnd))
    X = np.array([(p % (n + 1)) * h for p in free])
    Y = np.array([(p // (n + 1)) * h for p in free])
    return u[free], uex(X, Y), h, Araw


# ======================================================================
def errors(u, ue, h):
    e = u - ue
    return np.sqrt(np.sum(e**2) * h**2), np.max(np.abs(e))


def order(errs, hs):
    return [np.nan] + [np.log(errs[k - 1] / errs[k]) / np.log(hs[k - 1] / hs[k])
                       for k in range(1, len(errs))]


print("=" * 78)
print("E1.  The interior stencil, on a uniform mesh with constant k")
print("=" * 78)
n = 8
h = 1.0 / n
_, _, _, Afd = fd_solve(n)
_, _, _, Afv = fv_solve(n)
_, _, _, Ap1 = p1_solve(n)
_, _, _, Aq1 = q1_solve(n)

# representative interior rows, scaled so each approximates -Laplacian
def stencil_fd(A, n):
    m = n - 1
    i = j = m // 2
    p = (j - 1) * m + (i - 1)
    row = np.zeros((3, 3))
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            ii, jj = i + di, j + dj
            if 1 <= ii <= n - 1 and 1 <= jj <= n - 1:
                row[dj + 1, di + 1] = A[p, (jj - 1) * m + (ii - 1)]
    return row * h**2

def stencil_cell(A, n):
    i = j = n // 2
    p = j * n + i
    row = np.zeros((3, 3))
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            row[dj + 1, di + 1] = A[p, (j + dj) * n + (i + di)]
    return row          # already scaled: k*h/h, and RHS carries h^2

def stencil_node(A, n):
    i = j = n // 2
    nid = lambda a, b: b * (n + 1) + a
    p = nid(i, j)
    row = np.zeros((3, 3))
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            row[dj + 1, di + 1] = A[p, nid(i + di, j + dj)]
    return row

print("\nFinite difference, x h^2:\n", stencil_fd(Afd, n))
print("\nFinite volume (cell-centred), row as assembled:\n", stencil_cell(Afv, n))
print("\nFEM P1, right triangles, element stiffness row:\n", stencil_node(Ap1, n))
print("\nFEM Q1, bilinear quads, element stiffness row (x3 shown below):\n",
      stencil_node(Aq1, n))
print("\nFEM Q1 x 3:\n", 3 * stencil_node(Aq1, n))

print("\nmax|FD - FV| over the interior rows above :",
      np.max(np.abs(stencil_fd(Afd, n) - stencil_cell(Afv, n))))
print("max|FD - P1| over the interior rows above :",
      np.max(np.abs(stencil_fd(Afd, n) - stencil_node(Ap1, n))))
print("max|FD - Q1| over the interior rows above :",
      np.max(np.abs(stencil_fd(Afd, n) - stencil_node(Aq1, n))))

print()
print("=" * 78)
print("E2.  Convergence on the manufactured solution")
print("=" * 78)
ns = [8, 16, 32, 64, 128]
res = {}
for name, fn in (("FD  5-point", fd_solve), ("FVM cell-centred", fv_solve),
                 ("FEM P1 triangles", p1_solve), ("FEM Q1 quads", q1_solve)):
    L2, Li, hs, dofs = [], [], [], []
    for n in ns:
        u, ue, h, _ = fn(n)
        a, b = errors(u, ue, h)
        L2.append(a); Li.append(b); hs.append(h); dofs.append(len(u))
    res[name] = (hs, L2, Li, dofs)
    print(f"\n{name}")
    print(f"{'n':>5} {'dofs':>7} {'h':>10} {'L2 err':>12} {'ord':>6} {'Linf err':>12} {'ord':>6}")
    o2, oi = order(L2, hs), order(Li, hs)
    for k, n in enumerate(ns):
        print(f"{n:5d} {dofs[k]:7d} {hs[k]:10.5f} {L2[k]:12.4e} {o2[k]:6.2f} "
              f"{Li[k]:12.4e} {oi[k]:6.2f}")

import json
json.dump({k: {'h': v[0], 'L2': v[1], 'Linf': v[2], 'dofs': v[3]}
           for k, v in res.items()}, open(_data('conv.json'), 'w'))
