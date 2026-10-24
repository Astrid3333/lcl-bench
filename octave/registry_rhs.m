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
    case 'aizawa'
      a = params.a; b = params.b; c = params.c; d = params.d;
      e = params.e; fz = params.f;
      f = @(x) [(x(3) - b) * x(1) - d * x(2);
                d * x(1) + (x(3) - b) * x(2);
                c + a * x(3) - x(3).^3 / 3 ...
                - (x(1).^2 + x(2).^2) .* (1 + e * x(3)) ...
                + fz * x(3) * x(1).^3];
    case 'duffing'
      delta = params.delta; alpha = params.alpha; beta = params.beta;
      gamma = params.gamma; omega = params.omega;
      f = @(x) [x(2);
                -delta * x(2) - alpha * x(1) - beta * x(1)^3 ...
                + gamma * cos(x(3));
                omega];
    case 'fhn'
      a = params.a; b = params.b; delta = params.delta;
      r = params.r; omega0 = params.omega0;
      f = @(x) [(x(1) * (x(1) - a) * (1 - x(1)) - x(2)) / delta;
                x(1) - x(2) - b + r * sin(x(3));
                omega0];
    otherwise
      error('registry_rhs:unknown', 'Unknown ODE attractor: %s', name);
  end
end
