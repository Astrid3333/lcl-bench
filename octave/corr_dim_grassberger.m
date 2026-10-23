function dim = corr_dim_grassberger(pts, n_radii, c_lo, c_hi)
  % Grassberger-Procaccia correlation dimension of a point cloud.
  % C(r) = fraction of pairs closer than r; slope of log C vs log r in the
  % scaling window c_lo <= C <= c_hi. Pairwise Euclidean distances.
  if nargin < 2, n_radii = 30; end
  if nargin < 3, c_lo = 0.01; end
  if nargin < 4, c_hi = 0.5; end
  n = size(pts, 1);
  npairs = n * (n - 1) / 2;
  D = zeros(npairs, 1);
  k = 1;
  for i = 1:(n - 1)
    dif = pts(i + 1:end, :) - pts(i, :);
    dist = sqrt(sum(dif .^ 2, 2));
    D(k:(k + numel(dist) - 1)) = dist;
    k = k + numel(dist);
  end
  D = D(D > 0);
  if numel(D) < 10
    error('corr_dim_grassberger:pairs', 'too few nonzero pairs');
  end
  D = sort(D);
  r = logspace(log10(D(1)), log10(D(end)), n_radii);
  counts = zeros(n_radii, 1);
  for i = 1:n_radii
    counts(i) = sum(D <= r(i));
  end
  C = counts / numel(D);
  mask = C >= c_lo & C <= c_hi;
  if sum(mask) < 3
    error('corr_dim_grassberger:window', 'empty scaling window');
  end
  p = polyfit(log(r(mask)), log(C(mask)), 1);
  dim = p(1);
end
