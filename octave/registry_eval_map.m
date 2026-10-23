function res = registry_eval_map(e, outdir)
  % Evaluate a map entry: orbit reproducibility file + bit-exact head +
  % largest Lyapunov + optional Grassberger-Procaccia correlation dimension.
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
  has_ref = isfield(e, 'ref') && ~isempty(e.ref);

  % bit-exactness: first n iterates vs frozen reference sequence
  if has_ref && isfield(e.ref, 'orbit_head') && ~isempty(e.ref.orbit_head)
    oh = e.ref.orbit_head;
    ref_vals = oh.values;
    % jsondecode maps a nested float matrix to a plain double matrix
    head = orbit(1:oh.n, :);
    md = max(abs(head(:) - ref_vals(:)));
    tol = oh.abs_tol;
    checks(end + 1) = mkcheck('orbit_head_maxabs', md, 0, tol, md, md <= tol);
  end

  if has_ref && isfield(e.ref, 'lyap1') && ~isempty(e.ref.lyap1)
    val = e.ref.lyap1.value;
    tol = e.ref.lyap1.abs_tol;
    rel = abs(lam - val) / max(abs(val), eps);
    checks(end + 1) = mkcheck('lyap1', lam, val, tol, rel, abs(lam - val) <= tol);
  else
    checks(end + 1) = mkcheck('lyap1_observed', lam, NaN, NaN, NaN, true);
  end

  % Grassberger-Procaccia correlation dimension (protocol fixed in registry)
  if has_ref && isfield(e.ref, 'corr_dim') && ~isempty(e.ref.corr_dim)
    cd_ref = e.ref.corr_dim;
    proto = cd_ref.protocol;
    x = e.x0(:);
    for i = 1:proto.transient
      x = f(x);
    end
    pts = zeros(proto.n_points, e.dim);
    for i = 1:proto.n_points
      x = f(x);
      pts(i, :) = x';
    end
    dim = corr_dim_grassberger(pts, proto.n_radii, proto.c_lo, proto.c_hi);
    val = cd_ref.value;
    tol = cd_ref.rel_tol * abs(val);
    rel = abs(dim - val) / abs(val);
    checks(end + 1) = mkcheck('corr_dim', dim, val, tol, rel, rel <= cd_ref.rel_tol);
  end

  res = struct('name', e.name, 'family', e.family, 'pass', all([checks.pass]), ...
               'checks', checks, 'orbit_path', opath, 'lyap1_observed', lam);
end
