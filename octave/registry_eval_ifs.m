function res = registry_eval_ifs(e)
  % Evaluate an IFS entry: chaos-game box-counting dimension vs reference.
  pts = chaos_game(e);
  dim = box_counting_dim(pts, e.box_sizes(:)');
  checks = emptychecks();
  if isfield(e, 'ref') && isfield(e.ref, 'dim') && ~isempty(e.ref.dim)
    val = e.ref.dim.value;
    tol = e.ref.dim.rel_tol * abs(val);
    rel = abs(dim - val) / abs(val);
    checks(end + 1) = mkcheck('box_dim', dim, val, tol, rel, rel <= e.ref.dim.rel_tol);
  end
  if isempty(checks)
    checks = mkcheck('box_dim_observed', dim, NaN, NaN, NaN, true);
  end
  res = struct('name', e.name, 'family', e.family, 'pass', all([checks.pass]), ...
               'checks', checks, 'dim_observed', dim);
end

function pts = chaos_game(e)
  n = e.n_points;
  bi = e.burn_in;
  x = e.x0(:);
  pts = zeros(n, 2);
  if isfield(e.params, 'vertices')
    V = e.params.vertices;
    p = e.params.p(:)';
    cp = cumsum(p / sum(p));
    for i = 1:(n + bi)
      k = find(rand <= cp, 1);
      x = (x + V(k, :)') / 2;
      if i > bi
        pts(i - bi, :) = x';
      end
    end
  elseif isfield(e.params, 'maps')
    M = e.params.maps;
    pv = [M.p]';
    cp = cumsum(pv / sum(pv));
    for i = 1:(n + bi)
      k = find(rand <= cp, 1);
      m = M(k);
      x = [m.a * x(1) + m.b * x(2) + m.e;
           m.c * x(1) + m.d * x(2) + m.f];
      if i > bi
        pts(i - bi, :) = x';
      end
    end
  else
    error('registry_eval_ifs:params', 'Unknown IFS param shape for %s', e.name);
  end
end
