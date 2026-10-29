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


## Attractor diagnosis

`python/diagnose_attractor.py` compares a candidate implementation (a
degraded variant or an external trajectory) against the registry oracle
using the geometric fingerprint (λ1 + correlation dimension D2) and returns
`{diagnostico, potencial, confianza, evidence}`. Taxonomy: `equivalente`,
`colapso_por_redondeo`, `precision_insuficiente`, `artefacto_tipo_sna`,
`degradacion_de_orden`, `diferencia_por_integrador`,
`degradacion_de_precision`, `inestabilidad_numerica`,
`diferencia_no_explicada`. Thresholds are calibrated on measurements in
this repo (λ1 ulp band ±0.015 → LAM_NOISE 0.02; registry rel_tol 0.05 for
λ1; D2_REL 0.15 from the estimator's own noise: 12 consecutive oracle
windows give D2 1.663–1.828, std 3.3 %, max pairwise 9.4 % — a tight 8 %
false-mismatched 9 % of *same*-trajectory window pairs).

`--series` also measures λ1 from the trajectory alone
(`rosenstein_series`: z-scored metric so neutral ramp coordinates such as an
appended phase cannot dominate, Theiler window, median of per-pair
divergence slopes). Series λ1 is validated against the oracle only when the
neighbor distance is dense enough (`nn_med/D ≤ 5e-3`, the density gate) and
is compared under `LAM_REL_SERIES 0.10` (measured worst-case error 3.3 % on
realistic-size references); otherwise λ1 stays `None` with the reason in the
report and the verdict falls back to D2 under `D2_REL_SERIES 0.35`
(measured window max err 32.7 % on rossler at N=1000; ≤23.6 % at N=2000).
Scalar (1-column) series are
delay-embedded (m=3, τ=first ACF lag < 1/e) for λ1 only — D2 and
unique_frac stay on the raw rows. Use `--series-dt` to override the sample
spacing (default: 1.0 for maps, `dt_lyap` for ODEs).

```bash
python3 python/diagnose_attractor.py --system lorenz --variant euler
python3 python/diagnose_attractor.py --system henon --series some_orbit.txt
```

`python/parse_legacy_source.py` reads a legacy Fortran 77 / C source,
extracts language/precision/integrator/dt/system **with line evidence**, and
then **executes** it: `python/translate_f77.py` translates the strict F77
subset (fixed-form with labels/continuations, scalar and 1-D array
declarations with literal extents, arithmetic and intrinsics, block and
single-line `IF`, one `do` loop, `GOTO`/`GO TO` via a program-counter
dispatch, `write` after `end do`) and `python/translate_c.py` translates the
strict C subset into fixed-form F77, then reuses the same back-end (one
`for` loop with `<`/`<=` bounds, if/else-if/else chains, single-line `IF`,
`+=`/`++` increments, `!=`/`&&`/`||`/`!` operators, `f`-suffix literals;
refused constructs remap every error back to the original C line). Both run
under the precision declared in the source
(`REAL`/`float` → float32, `REAL*8`/`DOUBLE PRECISION`/`double` → float64),
so λ1 and D2 come
from the legacy recurrence itself (`candidate_source: executed_source`). The
translated program is stored in the report as `execution.translated_source`
and written next to `--out` as `*.translated.py`. Anything outside the
subset — several loops, non-literal array extents, labelled `DO`, other C
constructs — falls back to the
synthetic equivalent variant (`candidate_source: synthetic_equivalent`) with
an `execution.status` block (`executed_source` | `translation_failed` |
`unsupported_language`):

```bash
python3 python/parse_legacy_source.py fixtures/legacy/lorenz_euler_single.f
# single + euler + dt=0.01 -> executed_source -> precision_insuficiente (0.90)
python3 python/test_translate_c.py   # 127 branch checks for the C front-end
```

Fixture verdicts (all 6 executed, F77 and C): rk4/double →
`equivalente` (λ1 0.895 vs 0.904, D2 1.679 vs 1.827, both inside tolerance);
euler/single (F77 and C drivers, both λ1=nan) → `precision_insuficiente`; euler
dt=0.1 → `inestabilidad_numerica`; array-valued RK4 (`lorenz_rk4_array.f`,
structured path) → `equivalente`; block `IF` + `GOTO` euler
(`lorenz_euler_goto.f`, pc-dispatch path) → `degradacion_de_orden`
(λ1 1.047 vs 0.904, above the 0.0904 threshold). λ1 uses the legacy
program's own loop count as transient and the source's own `dt` for
normalisation.
Reports: `results/diagnose/*.json`, `results/legacy/*.json`
(+`*.translated.py`); `python3 python/check_verdicts.py` re-asserts the
preregistered verdict table (25 assertions / 13 artefacts).


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
| `python/diagnose_attractor.py` | fingerprint diagnosis (λ1+D2 vs oracle) |
| `python/parse_legacy_source.py` | legacy F77/C reader → evidence → verdict |
| `python/check_verdicts.py` | re-asserts the preregistered verdict table (25/25) |
| `python/translate_f77.py` | strict F77 subset → Python (executes legacy code) |
| `python/translate_c.py` | strict C subset → F77 → Python (executes legacy C) |
| `python/test_translate_c.py` | branch tests for the C front-end (127 checks) |
| `fixtures/legacy/` | 6 legacy-source fixtures (F77/C), incl. arrays + GOTO |
| `reports/informe-discrepancias.md` | discrepancy report (Spanish) |
| `run_all.sh` | one-command suite (5 steps, HDF5 optional) |
| `results/` | JSON results (timestamped), `discrepancy/`, `diagnose/`, `legacy/`, `registry_all.h5` |

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
