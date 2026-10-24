# Current State

- Goal: cross-language benchmark registry (chaotic attractors + physical potentials) with independent runners in **Octave, Python, and Fortran (gfortran)** that compute invariants and compare against cited references, plus a hub that diffs the three languages.
- Stage: **tasks 1–5 complete**. 18/18 entries green in all three languages, hub 18/18 ALL_PASS, discrepancy report + figures delivered (`reports/informe-discrepancias.md`), HDF5 export (`results/registry_all.h5`, schema registry_all/1), `./run_all.sh` = 5 steps (Octave → Python → Fortran → hub → HDF5, last one optional/skippable).
- Current task: cycle close (commit of tasks 4+5); then optional extensions.
- Next action: propose optional registry additions (Duffing, FitzHugh-Nagumo, De Jong/Lozi/Standard maps) + energy-drift integrator test (Verlet vs RK4) per design notes.
- Registries: `registries/attractors.json` (11 entries: logistic, henon, clifford, tinkerbell, sierpinski, barnsley_fern, lorenz, rossler, thomas, aizawa, hopfield_4), `registries/potentials.json` (7 entries: lennard_jones, morse, double_well, harmonic, hydrogen, yukawa, poschl_teller).
- Runners: `octave/registry_run.m`, `python/registry_run.py`, `fortran/frunner.f90` (program `f_runner`); selector semantics identical (all|attractors|potentials|family|entry-name). Hub: `python/hub_compare.py` (3-way: globs `_py_`, `_f_`, plain; writes `results/hub_compare.json`).
- Discrepancy tooling: `python/discrepancy_report.py` (36 quantities → `results/discrepancy/per_quantity.{json,csv}`), `python/sensitivity_aizawa.py` (14 scenarios), `python/figures_discrepancy.py` (fig1_spread + fig2_aizawa_sensitivity, openscience.mplstyle, venv matplotlib), `python/export_hdf5.py` (`.venv` + h5py).
- Key numbers: reference-backed spread median 1.78e-12, max 3.28e-3 (henon corr_dim); observed max 0.145 (aizawa lyap1, estimator-sensitive, ulp-proven); orbit-head residuals ≤2.22e-16 (bit-exact in 3 of 4 maps).
- Fortran: `kinds/fjson/fhandles/fkernels/frunner.f90` + `Makefile` (f_runner / f_runner_dbg / f_runner_asan); jcopy-based JSON deep copy; ASAN-clean 18/18.
- Tooling: gfortran 13.3.0; two gfortran bugs documented in `knowledge/fortran-gfortran-quirks.md`; `.venv/` (matplotlib 3.11.2, h5py 3.16.0) for figures/HDF5 — system python3 stays the canonical registry interpreter (numpy 1.26.4).
