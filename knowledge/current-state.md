# Current State

- Goal: cross-language benchmark registry (chaotic attractors + physical potentials) with a GNU Octave runner that computes invariants and compares against cited references.
- Stage: first deliverable complete — 16/16 entries PASS (`results/registry_run_20261023T022423.json`).
- Current task: cycle close (knowledge/notes/commit) after full-suite green run.
- Next action: user review; possible follow-ups are a Python runner hub (per original design) or adding entries.
- Registries: `registries/attractors.json` (10 entries: logistic, henon, clifford, tinkerbell, sierpinski, barnsley_fern, lorenz, rossler, thomas, hopfield_4), `registries/potentials.json` (6 entries: lennard_jones, morse, double_well, harmonic, hydrogen, yukawa).
- Runner: `octave/registry_run.m` with selector all|attractors|potentials|family|entry-name; JSON results in `results/`; orbit dumps in `results/orbits/`.
- Observed-only values (no citable ref yet): clifford lyap1≈0.2904, tinkerbell lyap1≈0.1995, thomas lyap1≈0.00647.
