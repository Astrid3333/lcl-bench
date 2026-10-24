# Cross-language benchmark registry

Chaotic attractors and physical potentials as JSON registry entries, evaluated
by **three independent implementations** (GNU Octave, Python/NumPy, Fortran/
gfortran) that each compare computed invariants against analytic or literature
references stored in the registry. A hub then diffs the three result sets.

## Quick start

```bash
./run_all.sh              # full suite + hub (builds Fortran automatically)
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
```

Results land in `results/registry_run_*.json` (Octave),
`results/registry_run_py_*.json` (Python),
`results/registry_run_f_*.json` (Fortran), and `results/hub_compare.json`.
Orbit dumps: `results/orbits/<name>_orbit.txt`.

## Layout

| Path | Role |
|---|---|
| `registries/attractors.json` | maps, IFS, ODEs, neural entries |
| `registries/potentials.json` | physical potential entries |
| `octave/` | Octave runner + kernels |
| `python/` | Python runner + kernels + hub |
| `fortran/` | Fortran runner + kernels + fjson parser + Makefile |
| `run_all.sh` | one-command suite (4 steps) |
| `results/` | JSON results (timestamped) |

## Invariants

- **Maps**: bit-exact first-20 iterates vs frozen `ref.orbit_head`
  (abs_tol 1e-12); largest Lyapunov (Benettin) vs `ref.lyap1`.
- **IFS**: box-counting dimension vs `ref.dim`.
- **ODEs**: largest Lyapunov (two-trajectory Benettin + RK4) vs `ref.lyap1`.
- **Neural (Hopfield)**: pattern fixed points.
- **Potentials**: analytic minima/forces, FD eigenvalues vs exact spectrum.

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
  ~15% (aizawa): a one-ulp change of `x0` alone moves the finite-time
  Lyapunov estimate by that much (see `results/hub_compare.json` field
  `obs_rel_diff`); reference-backed checks agree much tighter (≤4e-5
  rel). Pass flags are unaffected.
- Cosmetic: Python prints an absolute `Saved:` path, Octave/Fortran a
  relative one.
