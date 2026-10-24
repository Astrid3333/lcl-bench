function res = registry_eval_energy(e)
  % Energy-drift entry: velocity Verlet vs RK4 on a double-well potential,
  % comparing relative energy drift against frozen cross-language invariants.
  % Potential: V(x) = a*x^4 - b*x^2; F(x) = -dV/dx = -4*a*x^3 + 2*b*x.
  a = e.params.a;
  b = e.params.b;
  mass = e.mass;
  h = e.h;
  T = e.T;
  n = round(T / h);
  x = e.x0(1);
  v = e.v0(1);

  V = @(q) a * q.^4 - b * q.^2;
  F = @(q) -4.0 * a * q.^3 + 2.0 * b * q;

  E0 = 0.5 * mass * v * v + V(x);

  % velocity Verlet: kick-drift-kick
  xv = x; vv = v;
  max_rel = 0.0;
  for k = 1:n
    acc = F(xv) / mass;
    xv = xv + vv * h + 0.5 * acc * h * h;
    acc_new = F(xv) / mass;
    vv = vv + 0.5 * (acc + acc_new) * h;
    E = 0.5 * mass * vv * vv + V(xv);
    rel = abs(E - E0) / abs(E0);
    if rel > max_rel
      max_rel = rel;
    end
  end
  verlet_final = rel;
  verlet_max = max_rel;

  % RK4 on (x, v)
  f_y = @(y) [y(2); F(y(1)) / mass];
  y = [x; v];
  max_rel = 0.0;
  for k = 1:n
    k1 = f_y(y);
    k2 = f_y(y + h / 2.0 * k1);
    k3 = f_y(y + h / 2.0 * k2);
    k4 = f_y(y + h * k3);
    y = y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    E = 0.5 * mass * y(2)^2 + V(y(1));
    rel = abs(E - E0) / abs(E0);
    if rel > max_rel
      max_rel = rel;
    end
  end
  rk4_final = rel;
  rk4_max = max_rel;

  checks = emptychecks();
  ref = e.ref;
  keys = {'E0', 'verlet_final_rel_drift', 'verlet_max_rel_drift', ...
          'rk4_final_rel_drift', 'rk4_max_rel_drift'};
  vals = [E0, verlet_final, verlet_max, rk4_final, rk4_max];
  for j = 1:numel(keys)
    key = keys{j};
    computed = vals(j);
    if ~isfield(ref, key) || isempty(ref.(key))
      checks(end + 1) = mkcheck([key '_observed'], computed, NaN, NaN, NaN, true);
      continue;
    end
    r = ref.(key);
    val = r.value;
    if isfield(r, 'abs_tol')
      tol = r.abs_tol;
    else
      tol = r.rel_tol * abs(val);
    end
    rel = abs(computed - val) / max(abs(val), eps);
    checks(end + 1) = mkcheck(key, computed, val, tol, rel, ...
                              abs(computed - val) <= tol);
  end
  res = struct('name', e.name, 'family', e.family, ...
               'pass', all([checks.pass]), 'checks', checks);
end
