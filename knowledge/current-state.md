# Current State

- Goal: cross-language benchmark registry (chaotic attractors + physical potentials) with independent runners in Octave and Python that compute invariants and compare against cited references, plus a hub that diffs the two languages.
- Stage: 18/18 entries green in both languages; hub confirms 18/18 BOTH_PASS (`results/hub_compare.json`).
- Current task: task 3 (promote lessons to AGENTS.md); cycle close after commit.
- Next action: promote 1–2 verified lessons into AGENTS.md principles; then cycle close (notes + commit of new entries).
- Registries: `registries/attractors.json` (11 entries: logistic, henon, clifford, tinkerbell, sierpinski, barnsley_fern, lorenz, rossler, thomas, aizawa, hopfield_4), `registries/potentials.json` (7 entries: lennard_jones, morse, double_well, harmonic, hydrogen, yukawa, poschl_teller).
- Octave runner: `octave/registry_run.m` with selector all|attractors|potentials|family|entry-name; JSON results in `results/`; orbit dumps in `results/orbits/`.
- Python runner: `python/registry_run.py` (same selector semantics); hub: `python/hub_compare.py` diffs the latest Octave and Python JSON results.
- Observed-only values (no citable ref yet): clifford lyap1≈0.2904, tinkerbell lyap1≈0.1995, thomas lyap1≈0.00647, aizawa lyap1≈0.10733.
- Cross-language observed agreement (rel diff): aizawa 4e-15, thomas 2e-14, clifford 2e-16, tinkerbell 0. Same Benettin parameters give bit-close exponents; box-counting differs slightly (different RNG streams) but both within rel_tol.
