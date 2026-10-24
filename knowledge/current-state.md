# Current State

- Goal: cross-language benchmark registry (chaotic attractors + physical potentials) with independent runners in **Octave, Python, and Fortran (gfortran)** that compute invariants and compare against cited references, plus a hub that diffs the three languages.
- Stage: **18/18 entries green in all three languages**; hub 18/18 ALL_PASS `[PPP]`; `./run_all.sh` runs Octave → Python → build+Fortran → hub end-to-end (4 steps).
- Current task: cycle close (commit of Fortran runner + 3-language harness).
- Next action: discrepancy report + figures (task 4) — key evidence already gathered: observed-only spread (aizawa obs_rel_diff 0.145 traced to one-ulp x0 sensitivity of the finite-time Lyapunov estimator, d0-robust); optional HDF5 export (task 5); further entries (Duffing, FitzHugh-Nagumo, De Jong/Lozi/Standard), energy-drift integrator test.
- Registries: `registries/attractors.json` (11 entries: logistic, henon, clifford, tinkerbell, sierpinski, barnsley_fern, lorenz, rossler, thomas, aizawa, hopfield_4), `registries/potentials.json` (7 entries: lennard_jones, morse, double_well, harmonic, hydrogen, yukawa, poschl_teller).
- Runners: `octave/registry_run.m`, `python/registry_run.py`, `fortran/frunner.f90` (program `f_runner`); selector semantics identical (all|attractors|potentials|family|entry-name). Hub: `python/hub_compare.py` (3-way: globs `_py_`, `_f_`, plain; writes `results/hub_compare.json`).
- Fortran: `kinds/fjson/fhandles/fkernels/frunner.f90` + `Makefile` (f_runner / f_runner_dbg / f_runner_asan); jcopy-based JSON deep copy; ASAN-clean 18/18; console output structurally identical to Python (only `Saved:` path style differs).
- Invariants per family: maps = bit-exact first-20 iterates (`ref.orbit_head`, tol 1e-12) + Benettin lyap1 (+ Hénon GP corr_dim under frozen protocol); IFS = box-counting; ODE = Benettin lyap1 (RK4); neural = pattern fixed points; potentials = analytic minima/forces + FD eigenvalues.
- Observed-only values (no citable ref yet): clifford lyap1≈0.2904, tinkerbell lyap1≈0.1995, thomas lyap1≈0.00647, aizawa lyap1≈0.10733 (Octave; Fortran tinkerbell 0.1982, aizawa 0.0917 — finite-time estimator sensitivity, hub reports `obs_rel_diff`).
- Cross-language agreement: reference-backed checks rel ≤4e-5 (mostly ≤1e-9); observed: clifford 1.9e-16, thomas 2.1e-14, tinkerbell 6.4e-3, aizawa 1.45e-1; orbit_head maxdiff ≤2.2e-16.
- Tooling: gfortran 13.3.0; two gfortran bugs documented in `knowledge/fortran-gfortran-quirks.md`.
