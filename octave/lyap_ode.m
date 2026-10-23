function lambda = lyap_ode(f, x0, dt, t_transient, t_lyap, d0)
  % Largest Lyapunov exponent of a continuous system via two-trajectory
  % Benettin algorithm with fixed-step RK4 and per-step renormalization.
  % f: RHS handle; dt: renormalization step; t_transient/t_lyap: times.
  if nargin < 6, d0 = 1e-8; end
  h = dt;
  x = x0(:);
  n_tr = round(t_transient / h);
  for i = 1:n_tr
    x = rk4(x);
  end
  y = x + d0 * ones(size(x)) / norm(ones(size(x)));
  n = round(t_lyap / h);
  acc = 0;
  for i = 1:n
    x = rk4(x);
    y = rk4(y);
    d = y - x;
    nd = norm(d);
    acc = acc + log(nd / d0);
    y = x + d0 * (d / nd);
  end
  lambda = acc / (n * h);

  function xn = rk4(z)
    k1 = f(z);
    k2 = f(z + h / 2 * k1);
    k3 = f(z + h / 2 * k2);
    k4 = f(z + h * k3);
    xn = z + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
  end
end
