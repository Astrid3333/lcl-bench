# Cross-language benchmark registry

Chaotic attractors and physical potentials as JSON registry entries, evaluated
by **three independent implementations** (GNU Octave, Python/NumPy, Fortran/
gfortran) that each compare computed invariants against analytic or literature
references stored in the registry. A hub then diffs the three result sets.

## Quick start

```bash
./run_all.sh              # full suite + hub + HDF5 (builds Fortran automatically)
./run_all.sh map          # one family
./run_all.sh henon        # one entry
```

Or step by step:

```bash
octave-cli --no-gui --eval "addpath('octave'); registry_run('all');"
python3 python/registry_run.py all
make -s -C fortran f_runner
./fortran/f_runner all
python3 python/hub_compare.py
.venv/bin/python python/export_hdf5.py   # optional; needs h5py
```

Results land in `results/registry_run_*.json` (Octave),
`results/registry_run_py_*.json` (Python),
`results/registry_run_f_*.json` (Fortran), and `results/hub_compare.json`.
Orbit dumps: `results/orbits/<name>_orbit.txt`.
HDF5 bundle of everything: `results/registry_all.h5`.

## Discrepancy analysis

Where and why the three languages diverge is quantified in
`reports/informe-discrepancias.md` (with publication figures):

```bash
python3 python/discrepancy_report.py          # results/discrepancy/per_quantity.{json,csv}
python3 python/sensitivity_aizawa.py          # results/discrepancy/sensitivity_aizawa.json
.venv/bin/python python/figures_discrepancy.py
```

Headline: reference-backed checks agree to a median spread of 1.8e-12;
the largest observed gap (aizawa λ1, 15%) is finite-time-estimator
sensitivity to one-ulp x0 changes, demonstrated experimentally — not a bug.


## Layout

| Path | Role |
|---|---|
| `registries/attractors.json` | maps, IFS, ODEs, neural entries |
| `registries/potentials.json` | physical potential entries |
| `registries/energy_drift.json` | Verlet-vs-RK4 energy-drift test |
| `octave/` | Octave runner + kernels |
| `python/` | Python runner + kernels + hub |
| `fortran/` | Fortran runner + kernels + fjson parser + Makefile |
| `python/discrepancy_report.py` | per-quantity cross-language spread table |
| `python/sensitivity_aizawa.py` | aizawa λ1 sensitivity experiment |
| `python/figures_discrepancy.py` | publication figures for the report |
| `python/export_hdf5.py` | HDF5 bundle (needs h5py; `.venv` recommended) |
| `reports/informe-discrepancias.md` | discrepancy report (Spanish) |
| `run_all.sh` | one-command suite (5 steps, HDF5 optional) |
| `results/` | JSON results (timestamped), `discrepancy/`, `registry_all.h5` |

## Invariants

- **Maps**: bit-exact first-20 iterates vs frozen `ref.orbit_head`
  (abs_tol 1e-12); largest Lyapunov (Benettin) vs `ref.lyap1`.
- **IFS**: box-counting dimension vs `ref.dim`.
- **ODEs**: largest Lyapunov (two-trajectory Benettin + RK4) vs `ref.lyap1`.
- **Neural (Hopfield)**: pattern fixed points.
- **Potentials**: analytic minima/forces, FD eigenvalues vs exact spectrum.
- **Energy drift** (`double_well_verlet_rk4`): velocity-Verlet vs RK4 on
  V(x)=x⁴−x², x0=0.5, v0=0.4, h=0.01, T=1000; checks Verlet drift ≈0
  (rel_tol 1e-12), RK4 drift ≈0, and drift ratios (rel_tol 1e-5).

Hénon additionally reports the **Grassberger-Procaccia correlation
dimension** under the fixed protocol in `ref.corr_dim.protocol`.

Observed-only quantities (no citable reference yet) are reported with
`quantity*_observed` and pass=true; the hub compares them across languages.

## Adding an entry

1. Append a JSON object to the right registry (`family` ∈
   `map|ifs|ode|neural|potential`).
2. Provide `ref` fields with `value`, tolerance, and `source`
   (cite the literature; for exact formulas the source is the formula).
3. Run `./run_all.sh <name>` — no code changes needed unless the entry
   introduces a new equation (then add one branch to `octave/registry_rhs.m`
   or `registry_map_step.m`, the mirror in `python/handles.py`, and in
   `fortran/fhandles.f90`).

## Reproducibility notes

- All three languages compare against the **registry references**, not against
  each other; the hub only diffs the three result files afterwards.
- IFS chaos-game RNG differs per language (Python `default_rng(0)`,
  Octave unseeded, Fortran `random_seed` base 20261023+37j); box dims agree
  within rel_tol, not bit-exactly.
- Frozen map orbits in `ref.orbit_head` were generated once in float64 and
  cross-checked Octave↔Python at ≤1e-16 before freezing.
- Fortran build: `make -C fortran f_runner` (also `f_runner_dbg`,
  `f_runner_asan`, `clean`); flags include `-ffp-contract=off` so gfortran
  results stay IEEE-reproducible without FMA contraction.
- Observed-only values of chaotic systems show cross-language spread up to
  ~15% (aizawa) and ~11% (duffing): a one-ulp change of `x0` alone moves the
  finite-time Lyapunov estimate by that much (see `results/hub_compare.json`
  field `obs_rel_diff` and the sensitivity study in
  `reports/informe-discrepancias.md`); reference-backed checks agree much
  tighter (median 1.8e-12 rel). Pass flags are unaffected.
- Cosmetic: Python prints an absolute `Saved:` path, Octave/Fortran a
  relative one.
