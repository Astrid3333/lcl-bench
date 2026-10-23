function res = registry_eval_neural(e)
  % Evaluate a Hopfield entry: stored patterns must be fixed points of
  % sign(W p) with Hebbian zero-diagonal weights.
  P = e.params.patterns;
  n = e.params.n;
  W = (P' * P) / n;
  W = W - diag(diag(W));
  checks = emptychecks();
  for i = 1:size(P, 1)
    p = P(i, :)';
    s = sign(W * p);
    ok = all(s == p);
    checks(end + 1) = mkcheck(sprintf('fixed_point_pattern_%d', i), ...
                              double(ok), 1, 0, double(~ok), ok);
  end
  res = struct('name', e.name, 'family', e.family, 'pass', all([checks.pass]), ...
               'checks', checks);
end
