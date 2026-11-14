function res = registry_eval_cordic(e)
  % CORDIC (Volder 1959): rotacion de (1,0) por un angulo usando solo sumas,
  % restas y desplazamientos binarios, con tabla atan(2^-i) y ganancia K.
  % El vector roteado queda escalado por K: cos = x/K, sin = y/K.
  params = e.params;
  n_iter = double(params.iterations);
  mults = double(params.angle_multipliers);
  n_samples = double(params.samples);
  ref = e.ref;

  atan_tab = zeros(1, n_iter);
  K = 1.0;
  for i = 0:n_iter-1
    atan_tab(i+1) = atan(2^(-i));
    K = K * sqrt(1 + 2^(-2*i));
  end

  checks = emptychecks();

  K_ref = double(ref.gain_K.value);
  K_tol = double(ref.gain_K.abs_tol);
  checks(end+1) = mkcheck('gain_K', K, K_ref, K_tol, ...
                          abs(K - K_ref) / max(abs(K_ref), 1), ...
                          abs(K - K_ref) <= K_tol);

  cos_tol = double(ref.cos_err_max.abs_tol);
  sin_tol = double(ref.sin_err_max.abs_tol);
  res_tol = double(ref.residual_max.abs_tol);
  q15_tol = double(ref.q15_err_max.abs_tol);

  max_cos = 0.0;
  max_sin = 0.0;
  max_res = 0.0;
  max_q15 = 0.0;

  for k = 1:numel(mults)
    m = mults(k);
    theta = m * pi;
    [xr, yr, z] = cordic_rot(theta, n_iter, atan_tab);
    ec = abs(xr / K - cos(theta));
    es = abs(yr / K - sin(theta));
    tag = sprintf('m%g', m);
    checks(end+1) = mkcheck(['cos_err_' tag], ec, 0, cos_tol, ec, ec <= cos_tol);
    checks(end+1) = mkcheck(['sin_err_' tag], es, 0, sin_tol, es, es <= sin_tol);

    theta_f = fix(theta / pi * 32768);
    theta_r = (theta_f / 32768.0) * pi;
    [xq, yq] = cordic_rot(theta_r, n_iter, atan_tab);
    xq = fix(xq * 32768);
    yq = fix(yq * 32768);
    eq = max(abs(xq / (32768.0 * K) - cos(theta_r)), ...
             abs(yq / (32768.0 * K) - sin(theta_r)));
    max_q15 = max(max_q15, eq);
  end

  for k = 1:n_samples
    theta = -pi + 2 * pi * (k - 1) / (n_samples - 1);
    [xr, yr, z] = cordic_rot(theta, n_iter, atan_tab);
    max_cos = max(max_cos, abs(xr / K - cos(theta)));
    max_sin = max(max_sin, abs(yr / K - sin(theta)));
    max_res = max(max_res, abs(z));
  end

  checks(end+1) = mkcheck('cos_err_max', max_cos, 0, cos_tol, max_cos, ...
                          max_cos <= cos_tol);
  checks(end+1) = mkcheck('sin_err_max', max_sin, 0, sin_tol, max_sin, ...
                          max_sin <= sin_tol);
  checks(end+1) = mkcheck('residual_max', max_res, 0, res_tol, max_res, ...
                          max_res <= res_tol);
  checks(end+1) = mkcheck('q15_err_max', max_q15, 0, q15_tol, max_q15, ...
                          max_q15 <= q15_tol);

  res = struct('name', e.name, 'family', e.family, ...
               'pass', all([checks.pass]), 'checks', checks, ...
               'gain_K', K, 'cos_err_max', max_cos, 'sin_err_max', max_sin, ...
               'residual_max', max_res, 'q15_err_max', max_q15);
end

function [xr, yr, z] = cordic_rot(theta, n_iter, atan_tab)
  % Rotacion CORDIC: devuelve (K*cos(theta), K*sin(theta), residuo angular).
  x = 1.0;
  y = 0.0;
  t = theta;
  while t > pi
    t = t - 2 * pi;
  end
  while t < -pi
    t = t + 2 * pi;
  end
  sign_cos = 1.0;
  if t > pi / 2
    t = pi - t;
    sign_cos = -1.0;
  elseif t < -pi / 2
    t = -pi - t;
    sign_cos = -1.0;
  end
  z = t;
  for i = 1:n_iter
    sigma = 1.0;
    if z < 0
      sigma = -1.0;
    end
    f = 2^(-(i - 1));
    xn = x - sigma * y * f;
    yn = y + sigma * x * f;
    z = z - sigma * atan_tab(i);
    x = xn;
    y = yn;
  end
  xr = x * sign_cos;
  yr = y;
end
