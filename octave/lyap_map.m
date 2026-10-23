function lambda = lyap_map(f, x0, n, transient, d0)
  % Largest Lyapunov exponent of a discrete map via Benettin et al. (1980).
  % f: one-step handle; x0: column initial state; n: post-transient steps;
  % transient: discarded steps; d0: renormalization distance.
  if nargin < 5, d0 = 1e-8; end
  x = x0(:);
  for i = 1:transient
    x = f(x);
  end
  y = x + d0 * ones(size(x)) / norm(ones(size(x)));
  acc = 0;
  for i = 1:n
    x = f(x);
    y = f(y);
    d = y - x;
    nd = norm(d);
    acc = acc + log(nd / d0);
    y = x + d0 * (d / nd);
  end
  lambda = acc / n;
end
