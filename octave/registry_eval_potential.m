function res = registry_eval_potential(e)
  % Evaluate a potential entry against analytic references:
  % minima (LJ/Morse/double well), exact FD eigenvalues (harmonic/Coulomb),
  % analytic force vs central differences (Yukawa).
  checks = emptychecks();
  ref = e.ref;
  fopts = optimset('TolX', 1e-12);
  switch lower(e.name)
    case {'lennard_jones', 'morse', 'double_well'}
      switch lower(e.name)
        case 'lennard_jones'
          ep = e.params.eps; sg = e.params.sigma;
          V = @(r) 4 * ep * ((sg ./ r).^12 - (sg ./ r).^6);
          [xmin, vmin] = fminbnd(V, e.r_range(1), e.r_range(2), fopts);
        case 'morse'
          De = e.params.De; a = e.params.a; re = e.params.re;
          V = @(r) De * (1 - exp(-a * (r - re))).^2 - De;
          [xmin, vmin] = fminbnd(V, e.r_range(1), e.r_range(2), fopts);
        case 'double_well'
          a = e.params.a; b = e.params.b;
          V = @(x) a * x.^4 - b * x.^2;
          [xmin, vmin] = fminbnd(V, 0, e.r_range(2), fopts);  % positive well
      end
      if isfield(ref, 'r_min') && ~isempty(ref.r_min)
        val = ref.r_min.value; tol = ref.r_min.abs_tol;
        rel = abs(xmin - val) / abs(val);
        checks(end + 1) = mkcheck('r_min', xmin, val, tol, rel, abs(xmin - val) <= tol);
      end
      if isfield(ref, 'r_min_abs') && ~isempty(ref.r_min_abs)
        val = ref.r_min_abs.value; tol = ref.r_min_abs.abs_tol;
        rel = abs(abs(xmin) - val) / abs(val);
        checks(end + 1) = mkcheck('r_min_abs', abs(xmin), val, tol, rel, abs(abs(xmin) - val) <= tol);
      end
      if isfield(ref, 'V_min') && ~isempty(ref.V_min)
        val = ref.V_min.value; tol = ref.V_min.abs_tol;
        rel = abs(vmin - val) / max(abs(val), eps);
        checks(end + 1) = mkcheck('V_min', vmin, val, tol, rel, abs(vmin - val) <= tol);
      end
    case {'harmonic', 'hydrogen'}
      switch lower(e.name)
        case 'harmonic'
          g = linspace(e.x_grid.xmin, e.x_grid.xmax, e.x_grid.n)';
          Vg = 0.5 * e.params.k * g.^2;
          sigma = 0;  % spectrum [0.5, 1.5, ...]; shift 0 -> lowest states
        case 'hydrogen'
          % Half-line grid for u(r) = r R(r): Dirichlet at 0 and rmax
          % removes the cusp and the spurious near-origin FD states.
          n = e.x_grid.n;
          rmax = e.x_grid.xmax;
          h = rmax / (n + 1);
          g = h * (1:n)';
          Vg = -1 ./ g;
          sigma = -1;  % bound states negative; shift -1 -> ground state
      end
      k = numel(e.eigen_n);
      ev = fd_eigenvalues(g, Vg, k, sigma);
      for j = 1:k
        val = ref.eigenvalues.value(j);
        comp = ev(j);
        rel = abs(comp - val) / abs(val);
        checks(end + 1) = mkcheck(sprintf('E_%d', e.eigen_n(j)), comp, val, ...
                                  ref.eigenvalues.rel_tol * abs(val), rel, ...
                                  rel <= ref.eigenvalues.rel_tol);
      end
      if isfield(ref, 'omega') && ~isempty(ref.omega)
        w = sqrt(e.params.k / e.params.m);
        val = ref.omega.value; tol = ref.omega.abs_tol;
        checks(end + 1) = mkcheck('omega', w, val, tol, abs(w - val) / abs(val), abs(w - val) <= tol);
      end
    case 'yukawa'
      A = e.params.A; kappa = e.params.kappa;
      Vh = @(r) A * exp(-kappa * r) ./ r;
      Fh = @(r) A * exp(-kappa * r) .* (1 ./ r.^2 + kappa ./ r);
      hh = 1e-5;
      for r = e.r_check(:)'
        Fnum = -(Vh(r + hh) - Vh(r - hh)) / (2 * hh);
        Fan = Fh(r);
        val = Fan; tol = ref.force.abs_tol;
        rel = abs(Fnum - val) / abs(val);
        checks(end + 1) = mkcheck(sprintf('force_r_%g', r), Fnum, val, tol, rel, ...
                                  abs(Fnum - val) <= tol);
      end
      F1 = -(Vh(1 + hh) - Vh(1 - hh)) / (2 * hh);
      val = ref.force.value_at_r1; tol = ref.force.abs_tol;
      checks(end + 1) = mkcheck('force_r_1_ref', F1, val, tol, ...
                                abs(F1 - val) / abs(val), abs(F1 - val) <= tol);
    otherwise
      error('registry_eval_potential:name', 'Unknown potential: %s', e.name);
  end
  if isempty(checks)
    error('registry_eval_potential:nochecks', 'No checks produced for %s', e.name);
  end
  res = struct('name', e.name, 'family', e.family, 'pass', all([checks.pass]), ...
               'checks', checks);
end
