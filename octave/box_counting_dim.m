function dim = box_counting_dim(pts, sizes)
  % Box-counting dimension: slope of log N(eps) vs log(1/eps).
  % pts: n x d point cloud; sizes: vector of box edge lengths.
  lo = min(pts, [], 1);
  counts = zeros(numel(sizes), 1);
  for i = 1:numel(sizes)
    eps = sizes(i);
    q = floor((pts - lo) / eps);
    counts(i) = size(unique(q, 'rows'), 1);
  end
  good = counts > 1;
  p = polyfit(log(1 ./ sizes(good)), log(counts(good)), 1);
  dim = p(1);
end
