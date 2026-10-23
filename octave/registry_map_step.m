function f = registry_map_step(name, params)
  % Return a one-step handle f(x) (column vector) for a named discrete map.
  switch lower(name)
    case 'logistic'
      r = params.r;
      f = @(x) r * x(1) * (1 - x(1));
    case 'henon'
      a = params.a; b = params.b;
      f = @(x) [1 - a * x(1)^2 + x(2);
                b * x(1)];
    case 'clifford'
      a = params.a; b = params.b; c = params.c; d = params.d;
      f = @(x) [sin(a * x(2)) + c * cos(a * x(1));
                sin(b * x(1)) + d * cos(b * x(2))];
    case 'tinkerbell'
      a = params.a; b = params.b; c = params.c; d = params.d;
      f = @(x) [x(1)^2 - x(2)^2 + a * x(1) + b * x(2);
                2 * x(1) * x(2) + c * x(1) + d * x(2)];
    otherwise
      error('registry_map_step:unknown', 'Unknown map: %s', name);
  end
end
