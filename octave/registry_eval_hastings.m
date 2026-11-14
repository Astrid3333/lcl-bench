function res = registry_eval_hastings(e)
  % Polinomios de Hastings (1955) para sin, atan, exp y log1p.
  % Aproximaciones de Taylor de grado 5 con reduccion de argumento.
  params = e.params;
  n_samples = double(params.samples);
  sin_range = double(params.sin_range);
  atan_range = double(params.atan_range);
  exp_range = double(params.exp_range);
  log_range = double(params.log1p_range);
  ref = e.ref;

  checks = emptychecks();

  %% ---------- sin ----------
  err = range_max_err(@hastings_sin, @sin, sin_range(1), sin_range(2), n_samples);
  tol = double(ref.sin_err_max.abs_tol);
  checks(end+1) = mkcheck('sin_err_max', err, 0, tol, err, err <= tol);

  e6 = abs(hastings_sin(pi / 6) - sin(pi / 6));
  tol = double(ref.sin_err_pi6.abs_tol);
  checks(end+1) = mkcheck('sin_err_pi6', e6, 0, tol, e6, e6 <= tol);

  e2 = abs(hastings_sin(pi / 2) - sin(pi / 2));
  tol = double(ref.sin_err_pi2.abs_tol);
  checks(end+1) = mkcheck('sin_err_pi2', e2, 0, tol, e2, e2 <= tol);

  v = hastings_sin(pi / 2);
  vref = double(ref.sin_val_pi2.value);
  vtol = double(ref.sin_val_pi2.abs_tol);
  checks(end+1) = mkcheck('sin_val_pi2', v, vref, vtol, ...
                          abs(v - vref) / max(abs(vref), 1), abs(v - vref) <= vtol);

  %% ---------- atan ----------
  err = range_max_err(@hastings_atan, @atan, atan_range(1), atan_range(2), n_samples);
  tol = double(ref.atan_err_max.abs_tol);
  checks(end+1) = mkcheck('atan_err_max', err, 0, tol, err, err <= tol);

  e1 = abs(hastings_atan(1.0) - atan(1.0));
  tol = double(ref.atan_err_one.abs_tol);
  checks(end+1) = mkcheck('atan_err_one', e1, 0, tol, e1, e1 <= tol);

  %% ---------- exp ----------
  err = 0.0;
  for k = 1:n_samples
    x = exp_range(1) + (exp_range(2) - exp_range(1)) * (k - 1) / (n_samples - 1);
    err = max(err, abs(hastings_exp(x) - exp(x)) / exp(x));
  end
  tol = double(ref.exp_rel_err_max.abs_tol);
  checks(end+1) = mkcheck('exp_rel_err_max', err, 0, tol, err, err <= tol);

  e1 = abs(hastings_exp(1.0) - exp(1.0));
  tol = double(ref.exp_err_one.abs_tol);
  checks(end+1) = mkcheck('exp_err_one', e1, 0, tol, e1, e1 <= tol);

  v = hastings_exp(1.0);
  vref = double(ref.exp_val_one.value);
  vtol = double(ref.exp_val_one.abs_tol);
  checks(end+1) = mkcheck('exp_val_one', v, vref, vtol, ...
                          abs(v - vref) / max(abs(vref), 1), abs(v - vref) <= vtol);

  %% ---------- log1p ----------
  err = range_max_err(@hastings_log1p, @log1p, log_range(1), log_range(2), n_samples);
  tol = double(ref.log1p_err_max.abs_tol);
  checks(end+1) = mkcheck('log1p_err_max', err, 0, tol, err, err <= tol);

  e1 = abs(hastings_log1p(1.0) - log1p(1.0));
  tol = double(ref.log1p_err_one.abs_tol);
  checks(end+1) = mkcheck('log1p_err_one', e1, 0, tol, e1, e1 <= tol);

  v = hastings_log1p(1.0);
  vref = double(ref.log1p_val_one.value);
  vtol = double(ref.log1p_val_one.abs_tol);
  checks(end+1) = mkcheck('log1p_val_one', v, vref, vtol, ...
                          abs(v - vref) / max(abs(vref), 1), abs(v - vref) <= vtol);

  res = struct('name', e.name, 'family', e.family, ...
               'pass', all([checks.pass]), 'checks', checks);
end

function err = range_max_err(f, g, a, b, n)
  err = 0.0;
  for k = 1:n
    x = a + (b - a) * (k - 1) / (n - 1);
    err = max(err, abs(f(x) - g(x)));
  end
end

function y = hastings_poly(c, x)
  y = c(end);
  for i = length(c)-1:-1:1
    y = y * x + c(i);
  end
end

function y = hastings_sin(x)
  c = [1.0, -1/6, 1/120, -1/5040, 1/362880, -1/39916800];
  sign = 1.0;
  if x < 0
    sign = -1.0;
    x = -x;
  end
  k = floor(x / (2 * pi));
  x = x - k * 2 * pi;
  if x > pi
    x = pi - x;
  end
  if x > pi / 2
    x = pi - x;
  end
  y = sign * x * hastings_poly(c, x * x);
end

function y = hastings_atan(x)
  c = [1.0, -1/3, 1/5, -1/7, 1/9, -1/11];
  sign = 1.0;
  if x < 0
    sign = -1.0;
    x = -x;
  end
  if x >= 0.5 && x <= 2.0
    t = (x - 1.0) / (x + 1.0);
    y = sign * (pi / 4 + t * hastings_poly(c, t * t));
    return;
  end
  if x > 1.0
    inv = 1.0 / x;
    y = sign * (pi / 2 - inv * hastings_poly(c, inv * inv));
    return;
  end
  y = sign * x * hastings_poly(c, x * x);
end

function y = hastings_exp(x)
  c = [1.0, 1/2, 1/6, 1/24, 1/120, 1/720];
  if x > 1.0
    k = floor(x / 0.6931471805599453);
    rem = x - k * 0.6931471805599453;
    y = 2^k * (1.0 + rem * hastings_poly(c, rem));
    return;
  end
  if x < 0.0
    y = 1.0 / hastings_exp(-x);
    return;
  end
  k = 0;
  r = x;
  while r > 0.25
    r = r / 2.0;
    k = k + 1;
  end
  y = 1.0 + r * hastings_poly(c, r);
  for i = 1:k
    y = y * y;
  end
end

function y = hastings_log1p(x)
  c = [1.0, -1/2, 1/3, -1/4, 1/5, -1/6, 1/7, -1/8, 1/9, -1/10];
  if x > 1.0
    t = (x - 1.0) / (x + 1.0);
    y = 0.6931471805599453 + t * hastings_poly(c, t);
    return;
  end
  if x >= 0.25
    t = (x - 1.0) / 2.0;
    y = 0.6931471805599453 + t * hastings_poly(c, t);
    return;
  end
  y = x * hastings_poly(c, x);
end
