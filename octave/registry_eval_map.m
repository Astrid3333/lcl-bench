function res = registry_eval_map(e, outdir)
  % Evaluate a map entry: orbit reproducibility file + largest Lyapunov.
  here = fileparts(mfilename('fullpath'));
  if nargin < 2, outdir = fullfile(fileparts(here), 'results'); end
  f = registry_map_step(e.name, e.params);
  x = e.x0(:);
  norb = e.n_orbit;
  orbit = zeros(norb, e.dim);
  for i = 1:norb
    x = f(x);
    orbit(i, :) = x';
  end
  if ~all(isfinite(orbit(:)))
    error('registry_eval_map:nonfinite', 'Orbit of %s is non-finite', e.name);
  end
  odir = fullfile(outdir, 'orbits');
  if ~exist(odir, 'dir'), mkdir(odir); end
  opath = fullfile(odir, [e.name '_orbit.txt']);
  fid = fopen(opath, 'w');
  fmt = [repmat('%.17g ', 1, e.dim - 1), '%.17g\n'];
  for i = 1:norb
    fprintf(fid, fmt, orbit(i, :));
  end
  fclose(fid);

  lam = lyap_map(f, e.x0(:), e.n_lyap, e.transient);
  checks = emptychecks();
  if isfield(e, 'ref') && isfield(e.ref, 'lyap1') && ~isempty(e.ref.lyap1)
    val = e.ref.lyap1.value;
    tol = e.ref.lyap1.abs_tol;
    rel = abs(lam - val) / max(abs(val), eps);
    checks(end + 1) = mkcheck('lyap1', lam, val, tol, rel, abs(lam - val) <= tol);
  end
  if isempty(checks)
    checks = mkcheck('lyap1_observed', lam, NaN, NaN, NaN, true);
  end
  res = struct('name', e.name, 'family', e.family, 'pass', all([checks.pass]), ...
               'checks', checks, 'orbit_path', opath, 'lyap1_observed', lam);
end
