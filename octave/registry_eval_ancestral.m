function res = registry_eval_ancestral(e)
  % Algoritmos matematicos de tradiciones ancestrales: multiplicacion
  % egipcia, fracciones unitarias, raiz babilonica, seno de Aryabhata,
  % teorema chino del resto, sistema vigesimal maya y eliminacion china.
  params = e.params;
  mult = double(params.egyptian_multiply);
  frac = double(params.egyptian_fraction);
  sqrtp = params.babylonian_sqrt;
  aryp = params.aryabhata;
  tcrp = params.tcr;
  vig = double(params.vigesimal);
  lins = params.linear_system;
  ref = e.ref;

  checks = emptychecks();

  %% ---------- multiplicacion egipcia ----------
  prod = egyptian_mult(mult(1), mult(2));
  checks(end+1) = value_check('egyptian_product', prod, ref, 'egyptian_product');
  prod_neg = egyptian_mult(-mult(1), mult(2));
  checks(end+1) = value_check('egyptian_product_signed', prod_neg, ref, ...
                              'egyptian_product_signed');

  %% ---------- fracciones unitarias ----------
  [n_terms, dens] = egyptian_fraction(frac(1), frac(2), 32);
  for k = 1:min(n_terms, 2)
    checks(end+1) = value_check(sprintf('egyptian_den_%d', k), dens(k), ref, ...
                                sprintf('egyptian_den_%d', k));
  end
  vref = double(ref.egyptian_terms.value);
  tol = double(ref.egyptian_terms.abs_tol);
  checks(end+1) = mkcheck('egyptian_terms', n_terms, vref, tol, ...
                          abs(n_terms - vref) / max(abs(vref), 1), ...
                          abs(n_terms - vref) <= tol);

  s = 0.0;
  for k = 1:n_terms
    s = s + 1.0 / dens(k);
  end
  err = abs(s - frac(1) / frac(2));
  tol = double(ref.egyptian_sum_err.abs_tol);
  checks(end+1) = mkcheck('egyptian_sum_err', err, 0, tol, err, err <= tol);

  %% ---------- raiz cuadrada babilonica (Heron) ----------
  sq = babilonian_sqrt(double(sqrtp.x), double(sqrtp.iterations));
  err = abs(sq - sqrt(double(sqrtp.x)));
  tol = double(ref.babylonian_sqrt_err.abs_tol);
  checks(end+1) = mkcheck('babylonian_sqrt_err', err, 0, tol, err, err <= tol);

  %% ---------- seno de Aryabhata ----------
  theta = double(aryp.theta);
  pasos = double(aryp.steps);
  s_ary = aryabhata_sin(theta, pasos);
  paso = 2 * pi / pasos;
  indice = floor(theta / paso);
  err_tab = abs(s_ary - sin(indice * paso));
  tol = double(ref.aryabhata_table_err.abs_tol);
  checks(end+1) = mkcheck('aryabhata_table_err', err_tab, 0, tol, err_tab, ...
                          err_tab <= tol);
  err_sin = abs(s_ary - sin(theta));
  tol = double(ref.aryabhata_sin_err.abs_tol);
  checks(end+1) = mkcheck('aryabhata_sin_err', err_sin, 0, tol, err_sin, ...
                          err_sin <= tol);

  %% ---------- teorema chino del resto ----------
  residuos = double(tcrp.residues);
  modulos = double(tcrp.moduli);
  [x_crt, M] = chinese_remainder(residuos, modulos);
  checks(end+1) = value_check('tcr_x', x_crt, ref, 'tcr_x');
  checks(end+1) = value_check('tcr_M', M, ref, 'tcr_M');

  %% ---------- sistema vigesimal maya ----------
  [n_dig, digs] = vigesimal(vig, 16);
  for k = 1:min(n_dig, 3)
    checks(end+1) = value_check(sprintf('vigesimal_d%d', k - 1), digs(k), ref, ...
                                sprintf('vigesimal_d%d', k - 1));
  end
  vref = double(ref.vigesimal_digits.value);
  tol = double(ref.vigesimal_digits.abs_tol);
  checks(end+1) = mkcheck('vigesimal_digits', n_dig, vref, tol, ...
                          abs(n_dig - vref) / max(abs(vref), 1), ...
                          abs(n_dig - vref) <= tol);
  vrecon = 0.0;
  pot = 1.0;
  for k = 1:n_dig
    vrecon = vrecon + digs(k) * pot;
    pot = pot * 20;
  end
  checks(end+1) = value_check('vigesimal_value', vrecon, ref, 'vigesimal_value');

  %% ---------- eliminacion china (eliminacion gaussiana) ----------
  A = to_matrix(lins.A);
  b = to_vector(lins.b);
  n = numel(b);
  x = chinese_elimination(A, b);
  for k = 1:n
    checks(end+1) = value_check(sprintf('linear_x%d', k), x(k), ref, ...
                                sprintf('linear_x%d', k));
  end
  resid = 0.0;
  for i = 1:n
    acc = 0.0;
    for j = 1:n
      acc = acc + A(i, j) * x(j);
    end
    resid = max(resid, abs(acc - b(i)));
  end
  tol = double(ref.linear_residual_max.abs_tol);
  checks(end+1) = mkcheck('linear_residual_max', resid, 0, tol, resid, ...
                          resid <= tol);

  res = struct('name', e.name, 'family', e.family, ...
               'pass', all([checks.pass]), 'checks', checks);
end

function c = value_check(label, computed, ref, refname)
  vref = double(ref.(refname).value);
  tol = double(ref.(refname).abs_tol);
  c = mkcheck(label, computed, vref, tol, ...
              abs(computed - vref) / max(abs(vref), 1), ...
              abs(computed - vref) <= tol);
end

function M = to_matrix(v)
  % jsondecode puede devolver una matriz numerica o una celda de celdas.
  if iscell(v)
    M = zeros(rows(v), columns(v));
    for i = 1:rows(v)
      row = v{i};
      for j = 1:numel(row)
        M(i, j) = row{j};
      end
    end
  else
    M = double(v);
  end
end

function b = to_vector(v)
  if iscell(v)
    b = zeros(1, numel(v));
    for i = 1:numel(v)
      b(i) = v{i};
    end
  else
    b = double(v(:).');
  end
end

function y = egyptian_mult(a, b)
  res = 0.0;
  neg = false;
  if a < 0
    a = -a;
    neg = ~neg;
  end
  if b < 0
    b = -b;
    neg = ~neg;
  end
  while b > 0
    if mod(b, 2) == 1
      res = res + a;
    end
    a = a * 2;
    b = floor(b / 2);
  end
  if neg
    y = -res;
  else
    y = res;
  end
end

function [n, dens] = egyptian_fraction(num, den, max_terms)
  dens = zeros(1, max_terms);
  n = 0;
  while num > 0 && n < max_terms
    m = ceil(den / num);
    dens(n + 1) = m;
    n = n + 1;
    nn = num * m - den;
    nd = den * m;
    g = gcd(nn, nd);
    num = nn / g;
    den = nd / g;
  end
end

function [n, digs] = vigesimal(valor, max_digitos)
  % Convierte un entero a base 20 (digitos del menos al mas significativo).
  digs = zeros(1, max_digitos);
  n = 0;
  v = valor;
  if v == 0
    if max_digitos > 0
      digs(1) = 0;
      n = 1;
    end
    return;
  end
  while v > 0 && n < max_digitos
    n = n + 1;
    digs(n) = mod(v, 20);
    v = floor(v / 20);
  end
end

function y = babilonian_sqrt(x, iterations)
  if x < 0
    y = -1.0;
    return;
  end
  if x == 0
    y = 0.0;
    return;
  end
  v = x / 2.0;
  if v == 0
    v = 1.0;
  end
  for i = 1:iterations
    v = 0.5 * (v + x / v);
  end
  y = v;
end

function y = aryabhata_sin(theta, pasos)
  while theta < 0
    theta = theta + 2 * pi;
  end
  while theta > 2 * pi
    theta = theta - 2 * pi;
  end
  n = pasos;
  if n <= 0
    n = 360;
  end
  paso = 2 * pi / n;
  indice = floor(theta / paso);
  if indice < 0
    indice = 0;
  end
  if indice > n
    indice = n;
  end
  k = 2.0 * (1.0 - cos(paso));
  s_prev = 0.0;
  s_actual = sin(paso);
  if indice == 0
    y = 0.0;
    return;
  end
  if indice == 1
    y = s_actual;
    return;
  end
  for n = 1:indice-1
    s_siguiente = 2.0 * s_actual - s_prev - k * s_actual;
    s_prev = s_actual;
    s_actual = s_siguiente;
  end
  y = s_actual;
end

function [x, M] = chinese_remainder(residuos, modulos)
  n = numel(residuos);
  M = 1.0;
  for i = 1:n
    M = M * modulos(i);
  end
  x = 0.0;
  for i = 1:n
    Mi = M / modulos(i);
    inv = modular_inverse(mod(Mi, modulos(i)), modulos(i));
    x = x + residuos(i) * Mi * inv;
  end
  x = mod(x, M);
end

function inv = modular_inverse(a, m)
  [g, x] = extended_gcd(mod(a, m), m);
  if g ~= 1
    inv = -1;
    return;
  end
  inv = mod(x, m);
end

function [g, x, y] = extended_gcd(a, b)
  if b == 0
    g = a;
    x = 1;
    y = 0;
    return;
  end
  [g, x1, y1] = extended_gcd(b, mod(a, b));
  x = y1;
  y = x1 - floor(a / b) * y1;
end

function x = chinese_elimination(A, b)
  n = size(A, 1);
  A = A;
  b = b(:).';
  for col = 1:n
    piv = col;
    maxv = abs(A(col, col));
    for fila = col+1:n
      if abs(A(fila, col)) > maxv
        maxv = abs(A(fila, col));
        piv = fila;
      end
    end
    if maxv < 1e-14
      x = zeros(1, n);
      return;
    end
    if piv ~= col
      tmp = A(col, :);
      A(col, :) = A(piv, :);
      A(piv, :) = tmp;
      tmpb = b(col);
      b(col) = b(piv);
      b(piv) = tmpb;
    end
    for fila = col+1:n
      factor = A(fila, col) / A(col, col);
      for k = col:n
        A(fila, k) = A(fila, k) - factor * A(col, k);
      end
      b(fila) = b(fila) - factor * b(col);
    end
  end
  x = zeros(1, n);
  for i = n:-1:1
    acc = b(i);
    for j = i+1:n
      acc = acc - A(i, j) * x(j);
    end
    x(i) = acc / A(i, i);
  end
end
