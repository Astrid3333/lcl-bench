function res = registry_eval_lcg(e)
  % Calidad espectral de los LCG de Lehmer: secuencia exacta, planos de
  % Marsaglia en 3D, normales enteras minimas y periodo del ciclo.
  params = e.params;
  randu = params.randu;
  minstd = params.minstd;
  box_h = double(params.normal_box_h);
  mcoef = double(params.marsaglia_coeffs);
  n_samples = double(params.marsaglia_samples);
  pc = params.period_case;
  pcm = params.period_case_minstd;
  ref = e.ref;

  checks = emptychecks();

  %% ---------- secuencia exacta ----------
  seq_r = lcg_seq(double(randu.a), double(randu.c), double(randu.m), 1, 4);
  seq_m = lcg_seq(double(minstd.a), double(minstd.c), double(minstd.m), 1, 4);
  for k = 1:4
    checks(end+1) = value_check(sprintf('randu_x%d', k), seq_r(k), ...
                                ref, sprintf('randu_x%d', k));
    checks(end+1) = value_check(sprintf('minstd_x%d', k), seq_m(k), ...
                                ref, sprintf('minstd_x%d', k));
  end

  %% ---------- planos de Marsaglia (solo RANDU: 15 planos en 3D) ----------
  m_r = double(randu.m);
  seq = lcg_seq(double(randu.a), double(randu.c), m_r, 1, n_samples + 2);
  ks = zeros(1, n_samples);
  max_dev = 0.0;
  for i = 1:n_samples
    k = (mcoef(1) * seq(i) + mcoef(2) * seq(i + 1) + mcoef(3) * seq(i + 2)) / m_r;
    ks(i) = round(k);
    max_dev = max(max_dev, abs(k - ks(i)));
  end
  n_planes = numel(unique(ks));

  tol = double(ref.randu_plane_count.abs_tol);
  vref = double(ref.randu_plane_count.value);
  checks(end+1) = mkcheck('randu_plane_count', n_planes, vref, tol, ...
                          abs(n_planes - vref) / max(abs(vref), 1), ...
                          abs(n_planes - vref) <= tol);

  tol = double(ref.marsaglia_integrality.abs_tol);
  checks(end+1) = mkcheck('marsaglia_integrality', max_dev, 0, tol, max_dev, ...
                          max_dev <= tol);

  %% ---------- normales enteras minimas ----------
  [n_norm_r, min_norm_r] = lcg_normals(double(randu.a), m_r, box_h);
  [n_norm_m, ~] = lcg_normals(double(minstd.a), double(minstd.m), box_h);

  vref = double(ref.randu_normal_count.value);
  tol = double(ref.randu_normal_count.abs_tol);
  checks(end+1) = mkcheck('randu_normal_count', n_norm_r, vref, tol, ...
                          abs(n_norm_r - vref) / max(abs(vref), 1), ...
                          abs(n_norm_r - vref) <= tol);

  vref = double(ref.minstd_normal_count.value);
  tol = double(ref.minstd_normal_count.abs_tol);
  checks(end+1) = mkcheck('minstd_normal_count', n_norm_m, vref, tol, ...
                          abs(n_norm_m - vref) / max(abs(vref), 1), ...
                          abs(n_norm_m - vref) <= tol);

  vref = double(ref.randu_min_normal_norm.value);
  tol = double(ref.randu_min_normal_norm.abs_tol);
  checks(end+1) = mkcheck('randu_min_normal_norm', min_norm_r, vref, tol, ...
                          abs(min_norm_r - vref) / max(abs(vref), 1), ...
                          abs(min_norm_r - vref) <= tol);

  %% ---------- periodo del ciclo ----------
  per_r = lcg_period(double(pc.a), double(pc.c), double(pc.m), double(pc.x0));
  per_m = lcg_period(double(pcm.a), double(pcm.c), double(pcm.m), double(pcm.x0));

  vref = double(ref.randu_period.value);
  tol = double(ref.randu_period.abs_tol);
  checks(end+1) = mkcheck('randu_period', per_r, vref, tol, ...
                          abs(per_r - vref) / max(abs(vref), 1), ...
                          abs(per_r - vref) <= tol);

  vref = double(ref.minstd_period.value);
  tol = double(ref.minstd_period.abs_tol);
  checks(end+1) = mkcheck('minstd_period', per_m, vref, tol, ...
                          abs(per_m - vref) / max(abs(vref), 1), ...
                          abs(per_m - vref) <= tol);

  res = struct('name', e.name, 'family', e.family, ...
               'pass', all([checks.pass]), 'checks', checks, ...
               'randu_plane_count', n_planes, ...
               'randu_normal_count', n_norm_r, ...
               'minstd_normal_count', n_norm_m, ...
               'randu_min_normal_norm', min_norm_r, ...
               'randu_period', per_r, 'minstd_period', per_m);
end

function c = value_check(label, computed, ref, refname)
  vref = double(ref.(refname).value);
  tol = double(ref.(refname).abs_tol);
  c = mkcheck(label, computed, vref, tol, ...
              abs(computed - vref) / max(abs(vref), 1), ...
              abs(computed - vref) <= tol);
end

function seq = lcg_seq(a, c, m, x0, n)
  seq = zeros(1, n);
  x = mod(x0, m);
  for i = 1:n
    x = mod(a * x + c, m);
    seq(i) = x;
  end
end

function [count, min_norm] = lcg_normals(a, m, box)
  % Normales h = (h1,h2,h3) con h1 + h2*a + h3*a^2 = 0 (mod m) y |h_i| <= box.
  a2 = mod(a * a, m);
  count = 0;
  min_norm = Inf;
  for h3 = -box:box
    for h2 = -box:box
      h1 = mod(-(h2 * a + h3 * a2), m);
      if h1 > m / 2
        h1 = h1 - m;
      end
      if abs(h1) <= box && ~(h1 == 0 && h2 == 0 && h3 == 0)
        count = count + 1;
        nrm = sqrt(h1 * h1 + h2 * h2 + h3 * h3);
        if nrm < min_norm
          min_norm = nrm;
        end
      end
    end
  end
end

function p = lcg_period(a, c, m, x0)
  % Longitud del ciclo de x_{n+1} = (a*x_n + c) mod m.
  x = mod(x0, m);
  seen = false(1, m);
  p = 0;
  while ~seen(x + 1) && p < m
    seen(x + 1) = true;
    x = mod(a * x + c, m);
    p = p + 1;
  end
end
