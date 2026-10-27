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


def corr_dim_grassberger(pts, n_radii=30, c_lo=0.01, c_hi=0.5):
    """Grassberger-Procaccia correlation dimension of a point cloud.

    C(r) = fraction of pairs closer than r; slope of log C vs log r in the
    scaling window c_lo <= C <= c_hi. Pairwise Euclidean distances.
    """
    from scipy.spatial.distance import pdist
    pts = np.asarray(pts, dtype=float)
    d = pdist(pts)
    d = d[d > 0]
    if d.size < 10:
        raise RuntimeError("corr_dim: too few nonzero pairs")
    r = np.geomspace(d.min(), d.max(), n_radii)
    ds = np.sort(d)
    counts = np.searchsorted(ds, r, side="right")
    C = counts / float(ds.size)
    mask = (C >= c_lo) & (C <= c_hi)
    if int(mask.sum()) < 3:
        raise RuntimeError("corr_dim: empty scaling window")
    p = np.polyfit(np.log(r[mask]), np.log(C[mask]), 1)
    return float(p[0])


def acf_first_tau(x, max_lag=None):
    """First autocorrelation lag below 1/e (embedding delay choice)."""
    x = np.asarray(x, dtype=float).ravel()
    x = x - x.mean()
    n = x.size
    denom = float(np.dot(x, x))
    if denom <= 0.0:
        return 1
    max_lag = int(max_lag) if max_lag else max(1, min(n // 4, 500))
    for k in range(1, max_lag + 1):
        if np.dot(x[:-k], x[k:]) / denom < 1.0 / np.e:
            return k
    return max(1, max_lag)


def delay_embed(x, m=3, tau=None):
    """Delay-coordinate embedding: returns (rows, tau)."""
    x = np.asarray(x, dtype=float).ravel()
    tau = int(tau) if tau else acf_first_tau(x)
    tau = max(1, tau)
    rows = x.size - (m - 1) * tau
    if rows < 10:
        return np.empty((0, m)), tau
    cols = [x[i:i + rows] for i in range(0, m * tau, tau)]
    return np.column_stack(cols), tau


def lyap1_series(pts, h, window=20, sat_frac=0.2, t_min=0, min_points=6,
                 min_span=8, min_pairs=100, max_anchors=2000, t_max=5000,
                 density_gate=5e-3):
    """Largest Lyapunov exponent from a single dense trajectory (Rosenstein-style).

    Neighbor pairs (with a Theiler window) are followed in time; per-pair
    OLS slopes of ln d(t) up to the saturation diameter are aggregated by
    their median and divided by the sample spacing h. Columns are
    z-scored so neutral ramp coordinates (e.g. an appended phase) cannot
    dominate the metric. Returns (value or None, info dict); None means
    not measurable from this series, with info["reason"] saying why.
    """
    if not h or h <= 0:
        return None, {"reason": f"invalid sample spacing h={h}"}
    pts = np.asarray(pts, dtype=float)
    n = pts.shape[0]
    mu = np.nanmean(pts, axis=0)
    sd = np.nanstd(pts, axis=0)
    z = (pts - mu) / np.where(sd > 0, sd, 1.0)
    lo, hi = np.nanmin(z, axis=0), np.nanmax(z, axis=0)
    D = float(np.linalg.norm(hi - lo))
    if not np.isfinite(D) or D <= 0:
        return None, {"reason": "constant series (zero extent)"}
    T = int(min(t_max, n))
    stride = max(1, n // int(max_anchors))
    anchors = np.arange(0, n, stride)
    j = np.full(anchors.size, -1, dtype=np.int64)
    d0 = np.full(anchors.size, np.inf)
    for s in range(0, anchors.size, 256):
        a = anchors[s:s + 256]
        diff = z[a][:, None, :] - z[None, :, :]
        d = np.sqrt((diff * diff).sum(-1))
        d[np.abs(a[:, None] - np.arange(n)[None, :]) < window] = np.inf
        k = np.argmin(d, axis=1)
        dm = d[np.arange(a.size), k]
        j[s:s + a.size] = np.where(np.isfinite(dm), k, -1)
        d0[s:s + a.size] = dm
    ok = (j >= 0) & (d0 > 0)
    anchors, j, d0 = anchors[ok], j[ok], d0[ok]
    if anchors.size < min_pairs:
        return None, {"reason": f"too few neighbor pairs ({anchors.size} < {min_pairs})"}
    d0_med = float(np.median(d0))
    if d0_med / D > density_gate:
        return None, {"reason": f"density gate: nn_med/D={d0_med / D:.4f} > {density_gate}",
                      "d0D": d0_med / D}
    t = np.arange(T)
    dmat = np.full((anchors.size, T), np.nan)
    for s in range(0, anchors.size, 384):
        ii = anchors[s:s + 384]
        jj = j[s:s + 384]
        valid = (ii[:, None] + t[None, :] < n) & (jj[:, None] + t[None, :] < n)
        ic = np.clip(ii[:, None] + t[None, :], 0, n - 1)
        jc = np.clip(jj[:, None] + t[None, :], 0, n - 1)
        dmat[s:s + ii.size] = np.where(
            valid, np.sqrt(((z[ic] - z[jc]) ** 2).sum(-1)), np.nan)
    hit = (dmat >= sat_frac * D) | np.isnan(dmat)
    first = np.where(hit.any(1), np.argmax(hit, axis=1), T)
    M = (np.arange(T)[None, :] < first[:, None]) & (t[None, :] >= t_min) & \
        (dmat >= 1e-12 * D)
    Mf = M.astype(np.float64)
    y = np.where(M, np.log(np.where(M, dmat, 1.0)), 0.0)
    tt = t.astype(float)
    cnt = Mf.sum(1)
    sx, sxx = Mf @ tt, Mf @ (tt * tt)
    sy, sxy = y.sum(1), (y * tt).sum(1)
    idx = np.where(M, np.arange(T)[None, :], -1)
    firstv = np.where(cnt > 0, np.where(M, np.arange(T)[None, :], T).min(1), -1)
    span = idx.max(1) - firstv
    nk = cnt.astype(float)
    den = nk * sxx - sx * sx
    keep = (cnt >= min_points) & (span >= min_span) & (np.abs(den) > 1e-12)
    sl = (nk * sxy - sx * sy)[keep] / den[keep]
    sl = sl[np.isfinite(sl)]
    if sl.size < min_pairs:
        return None, {"reason": f"too few usable fits ({sl.size} < {min_pairs})",
                      "d0D": d0_med / D}
    lam = float(np.median(sl)) / h
    if not np.isfinite(lam):
        return None, {"reason": "non-finite slope estimate", "d0D": d0_med / D}
    return lam, {"d0D": d0_med / D, "npairs": int(sl.size)}


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
