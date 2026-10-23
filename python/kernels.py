"""Numerical kernels: Benettin Lyapunov, box-counting, FD eigenvalues, fixed-step RK4."""
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh


def lyap_map(f, x0, n, transient, d0=1e-8):
    """Largest Lyapunov exponent of a discrete map (Benettin et al. 1980)."""
    x = np.asarray(x0, dtype=float).copy()
    for _ in range(transient):
        x = f(x)
    ones = np.ones_like(x)
    y = x + d0 * ones / np.linalg.norm(ones)
    acc = 0.0
    for _ in range(n):
        x = f(x)
        y = f(y)
        d = y - x
        nd = np.linalg.norm(d)
        acc += np.log(nd / d0)
        y = x + d0 * (d / nd)
    return acc / n


def lyap_ode(f, x0, dt, t_transient, t_lyap, d0=1e-8):
    """Largest Lyapunov exponent via two-trajectory Benettin with fixed-step RK4."""
    h = dt

    def rk4(z):
        k1 = f(z)
        k2 = f(z + h / 2.0 * k1)
        k3 = f(z + h / 2.0 * k2)
        k4 = f(z + h * k3)
        return z + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    x = np.asarray(x0, dtype=float).copy()
    n_tr = int(round(t_transient / h))
    for _ in range(n_tr):
        x = rk4(x)
    ones = np.ones_like(x)
    y = x + d0 * ones / np.linalg.norm(ones)
    n = int(round(t_lyap / h))
    acc = 0.0
    for _ in range(n):
        x = rk4(x)
        y = rk4(y)
        d = y - x
        nd = np.linalg.norm(d)
        acc += np.log(nd / d0)
        y = x + d0 * (d / nd)
    return acc / (n * h)


def box_counting_dim(pts, sizes):
    """Box-counting dimension: slope of log N(eps) vs log(1/eps)."""
    pts = np.asarray(pts, dtype=float)
    sizes = np.asarray(sizes, dtype=float)
    lo = pts.min(axis=0)
    counts = np.zeros(sizes.size)
    for i, eps in enumerate(sizes):
        q = np.floor((pts - lo) / eps).astype(np.int64)
        counts[i] = np.unique(q, axis=0).shape[0]
    good = counts > 1
    logx = np.log(1.0 / sizes[good])
    logy = np.log(counts[good])
    # least-squares slope (same as polyfit deg 1)
    A = np.vstack([logx, np.ones_like(logx)]).T
    p = np.linalg.lstsq(A, logy, rcond=None)[0]
    return float(p[0])


def fd_eigenvalues(g, V, k, sigma):
    """Lowest k eigenvalues of -1/2 d2/dx2 + V on grid g via sparse FD (Dirichlet)."""
    g = np.asarray(g, dtype=float)
    V = np.asarray(V, dtype=float)
    h = g[1] - g[0]
    n = g.size
    c = 0.5 / h ** 2
    main = 2.0 * c * np.ones(n) + V
    off = -c * np.ones(n - 1)
    H = diags([off, main, off], [-1, 0, 1], format="csc")
    # eigsh with shift-invert: find eigenvalues near sigma
    ev = eigsh(H, k=k, sigma=sigma, which="LM", return_eigenvectors=False)
    return np.sort(ev)
