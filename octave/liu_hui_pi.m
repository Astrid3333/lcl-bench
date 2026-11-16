function out = liu_hui_pi(iterations, levels, stable)
  % Pi de Arquimedes (poligono circunscrito, hexagono inicial) + extrapolacion de Liu Hui/Romberg.
  % stable=true usa t/(sqrt(1+t^2)+1); false usa (sqrt(1+t^2)-1)/t (pierde digitos y colapsa).
  if nargin < 1, iterations = 10; end
  if nargin < 2, levels = 3; end
  if nargin < 3, stable = true; end
  t = 1/sqrt(3); n = 6; A = zeros(1, iterations+1);
  for k = 0:iterations
    A(k+1) = n*t;
    if stable, t = t/(sqrt(1+t^2)+1); else, t = (sqrt(1+t^2)-1)/t; end
    n = 2*n;
  end
  R = A;
  for j = 1:levels
    R = R(2:end) + (R(2:end)-R(1:end-1))/(4^j - 1);
  end
  out.pi = R(end); out.error = abs(R(end)-pi); out.sides = 6*2^iterations;
  out.archimedes_plain = A(end); out.levels = levels; out.stable = stable;
end
