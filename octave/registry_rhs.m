function f = registry_rhs(name, params)
  % Return an RHS handle f(x) (column vector) for a named ODE attractor.
  % Resolving the handle once avoids a switch dispatch on every RK4 eval.
  switch lower(name)
    case 'lorenz'
      s = params.sigma; r = params.rho; b = params.beta;
      f = @(x) [s * (x(2) - x(1));
                x(1) * (r - x(3)) - x(2);
                x(1) * x(2) - b * x(3)];
    case 'rossler'
      a = params.a; b = params.b; c = params.c;
      f = @(x) [-x(2) - x(3);
                x(1) + a * x(2);
                b + x(3) * (x(1) - c)];
    case 'thomas'
      b = params.b;
      f = @(x) [sin(x(2)) - b * x(1);
                sin(x(3)) - b * x(2);
                sin(x(1)) - b * x(3)];
    otherwise
      error('registry_rhs:unknown', 'Unknown ODE attractor: %s', name);
  end
end
