# gfortran quirks (verified in this workspace)

## 1. Recursive intrinsic assignment of derived types is SHALLOW (gfortran 13.3)

Assigning a derived type with allocatable components inside a `recursive`
procedure (`d = s`) copies pointers, not data: both variables free the same
buffers at frame exit → double free / use-after-free (ASAN). Plain procedures
deep-copy correctly. Fix: an explicit `jcopy(d, s)` recursive helper that
copies every field/component by component (`fjson.f90`, 39 call sites in
`frunner.f90`). Minimal repro: `/tmp/opencode/frepro/v1.f90` (fails) vs
`v4.f90` (passes).

## 2. `.and.` does not short-circuit

`do while (len_trim(s) > 0 .and. s(len_trim(s):len_trim(s)) == '0')` crashes
with "Substring out of bounds" when `len_trim(s)==0`: BOTH operands are
evaluated. Rewrite as an outer `do while` plus an early `exit`. Hit in the
`fg` %g formatter (`frunner.f90`, three strip loops).

## 3. IEEE determinism

Build with `-ffp-contract=off` (default `-O2` otherwise contracts to FMA);
then gfortran matches Python/NumPy op-for-op on identical source expressions
to ~1e-14 (clifford obs_rel_diff 1.9e-16). Chaotic finite-time estimates can
still diverge by ~15% from one-ulp input differences alone — that is
estimator sensitivity, not a codegen bug (see aizawa in README).
