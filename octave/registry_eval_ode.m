function res = registry_eval_ode(e)
  % Evaluate an ODE entry: largest Lyapunov exponent vs reference.
  f = registry_rhs(e.name, e.params);
  lam = lyap_ode(f, e.x0(:), e.dt_lyap, e.t_transient, e.t_lyap);
  if ~isfinite(lam)
    error('registry_eval_ode:nonfinite', 'Lyapunov of %s is non-finite', e.name);
  end
  checks = emptychecks();
  if isfield(e, 'ref') && isfield(e.ref, 'lyap1') && ~isempty(e.ref.lyap1)
    val = e.ref.lyap1.value;
    if isfield(e.ref.lyap1, 'abs_tol') && ~isempty(e.ref.lyap1.abs_tol)
      tol = e.ref.lyap1.abs_tol;
    else
      tol = e.ref.lyap1.rel_tol * abs(val);
    end
    rel = abs(lam - val) / max(abs(val), eps);
    checks(end + 1) = mkcheck('lyap1', lam, val, tol, rel, abs(lam - val) <= tol);
  end
  if isempty(checks)
    checks = mkcheck('lyap1_observed', lam, NaN, NaN, NaN, true);
  end
  res = struct('name', e.name, 'family', e.family, 'pass', all([checks.pass]), ...
               'checks', checks, 'lyap1_observed', lam);
end
