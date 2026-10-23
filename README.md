# Cross-language benchmark registry

Chaotic attractors and physical potentials as JSON registry entries, evaluated
by **two independent implementations** (GNU Octave and Python/NumPy) that each
compare computed invariants against analytic or literature references stored in
the registry. A hub then diffs the two result sets.

## Quick start

```bash
./run_all.sh              # full suite + hub
./run_all.sh map          # one family
./run_all.sh henon        # one entry
```

Or step by step:

```bash
octave-cli --no-gui --eval "addpath('octave'); registry_run('all');"
python3 python/registry_run.py all
python3 python/hub_compare.py
```

Results land in `results/registry_run_*.json` (Octave),
`results/registry_run_py_*.json` (Python), and `results/hub_compare.json`.
Orbit dumps: `results/orbits/<name>_orbit.txt`.

## Layout

| Path | Role |
|---|---|
| `registries/attractors.json` | maps, IFS, ODEs, neural entries |
| `registries/potentials.json` | physical potential entries |
| `octave/` | Octave runner + kernels |
| `python/` | Python runner + kernels + hub |
| `run_all.sh` | one-command suite |
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
   or `registry_map_step.m` and the mirror in `python/handles.py`).

## Reproducibility notes

- Both languages compare against the **registry references**, not against each
  other; the hub only diffs the two result files afterwards.
- IFS chaos-game uses a fixed RNG seed (Python: `default_rng(0)`; Octave:
  `rand('state')` is not seeded to match — box dims agree within rel_tol, not
  bit-exactly).
- Frozen map orbits in `ref.orbit_head` were generated once in float64 and
  cross-checked Octave↔Python at ≤1e-16 before freezing.
