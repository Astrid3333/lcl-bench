function ev = fd_eigenvalues(g, V, k, sigma)
  % Lowest k eigenvalues of -1/2 d2/dx2 + V(x) on grid g via sparse FD.
  % Dirichlet BCs at the grid ends. sigma: shift for eigs shift-invert.
  h = g(2) - g(1);
  n = numel(g);
  c = 0.5 / h^2;
  main = 2 * c * ones(n, 1) + V(:);
  off = -c * ones(n, 1);
  H = spdiags([off main off], -1:1, n, n);
  opts = struct('issym', true, 'isreal', true, 'tol', 1e-12, 'maxit', 5000);
  ev = sort(eigs(H, k, sigma, opts));
end
