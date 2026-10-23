# Current State

- Goal: cross-language benchmark registry (chaotic attractors + physical potentials) with independent runners in Octave and Python that compute invariants and compare against cited references, plus a hub that diffs the two languages.
- Stage: 18/18 entries green in both languages; hub 18/18 BOTH_PASS; one-command `./run_all.sh` runs Octave → Python → hub end-to-end.
- Current task: cycle close (commit of reproducibility + invariant work).
- Next action: remaining design gaps — more entries (Duffing, FitzHugh-Nagumo, De Jong/Lozi/Standard maps), energy-drift integrator test, a third language (gfortran available), discrepancy report + figures, HDF5 export (optional).
- Registries: `registries/attractors.json` (11 entries: logistic, henon, clifford, tinkerbell, sierpinski, barnsley_fern, lorenz, rossler, thomas, aizawa, hopfield_4), `registries/potentials.json` (7 entries: lennard_jones, morse, double_well, harmonic, hydrogen, yukawa, poschl_teller).
- Octave runner: `octave/registry_run.m` with selector all|attractors|potentials|family|entry-name; JSON results in `results/`; orbit dumps in `results/orbits/`.
- Python runner: `python/registry_run.py` (same selector semantics); hub: `python/hub_compare.py` diffs the latest Octave and Python JSON results.
- Invariants per family: maps = bit-exact first-20 iterates (`ref.orbit_head`, tol 1e-12) + Benettin lyap1 (+ Hénon GP corr_dim under frozen protocol); IFS = box-counting; ODE = Benettin lyap1 (RK4); neural = pattern fixed points; potentials = analytic minima/forces + FD eigenvalues.
- Observed-only values (no citable ref yet): clifford lyap1≈0.2904, tinkerbell lyap1≈0.1995, thomas lyap1≈0.00647, aizawa lyap1≈0.10733.
- Cross-language agreement (rel diff): aizawa 4e-15, thomas 2e-14, clifford 2e-16, tinkerbell 0; orbit_head maxdiff ≤2.2e-16; corr_dim (henon) identical to 12 digits (1.20780802419).
- Tooling: `gfortran 13.3.0` available for a future third-language runner; `run_all.sh` is the single entry point.
