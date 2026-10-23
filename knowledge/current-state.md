# Current State

- Goal: cross-language benchmark registry (chaotic attractors + physical potentials) with independent runners in Octave and Python that compute invariants and compare against cited references, plus a hub that diffs the two languages.
- Stage: Python hub complete — both languages independently PASS 16/16 on the same registry; hub confirms 16/16 BOTH_PASS (`results/hub_compare.json`).
- Current task: cycle close (knowledge/notes/commit) after green hub run.
- Next action: add more registry entries (task 2), then promote lessons to AGENTS.md principles (task 3).
- Registries: `registries/attractors.json` (10 entries: logistic, henon, clifford, tinkerbell, sierpinski, barnsley_fern, lorenz, rossler, thomas, hopfield_4), `registries/potentials.json` (6 entries: lennard_jones, morse, double_well, harmonic, hydrogen, yukawa).
- Octave runner: `octave/registry_run.m` with selector all|attractors|potentials|family|entry-name; JSON results in `results/`; orbit dumps in `results/orbits/`.
- Python runner: `python/registry_run.py` (same selector semantics); hub: `python/hub_compare.py` diffs the latest Octave and Python JSON results.
- Observed-only values (no citable ref yet): clifford lyap1≈0.2904, tinkerbell lyap1≈0.1995, thomas lyap1≈0.00647.
- Cross-language observed agreement (rel diff): thomas 2e-14, clifford 2e-16, tinkerbell 0. Same Benettin parameters give bit-close exponents; box-counting differs slightly (different RNG streams) but both within rel_tol.
