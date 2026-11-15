function res = registry_eval_prng_legacy(e)
  % PRNG obsoletos: medio del cuadrado de von Neumann (1946) y
  % Wichmann-Hill (Algorithm AS 183, 1982).
  params = e.params;
  ref = e.ref;

  checks = emptychecks();

  if strcmp(e.name, 'middle_square_vn')
    digits = double(params.digits);
    seed = double(params.seed);
    nseq = double(params.seq_values);
    n_states = 10 ^ digits;
    x = seed;
    for k = 1:nseq
      x = ms_step(x, digits);
      checks(end+1) = value_check(sprintf('seq_%d', k), x, ref, ...
                                  sprintf('seq_%d', k));
    end

    fixed = [];
    for s = 0:(n_states - 1)
      if ms_step(s, digits) == s
        fixed(end+1) = s;
      end
    end
    nonzero = fixed(fixed != 0);
    checks(end+1) = value_check('fixed_points_nonzero_count', numel(nonzero), ...
                                ref, 'fixed_points_nonzero_count');
    checks(end+1) = value_check('zero_is_fixed', any(fixed == 0), ref, ...
                                'zero_is_fixed');

    outs = zeros(1, n_states);
    for s = 0:(n_states - 1)
      outs(s + 1) = ms_step(s, digits);
    end
    checks(end+1) = value_check(sprintf('distinct_outputs_n%d', digits), ...
                                numel(unique(outs)), ref, ...
                                sprintf('distinct_outputs_n%d', digits));

    seen = false(1, n_states);
    x = seed;
    p = 0;
    while ~seen(x + 1)
      seen(x + 1) = true;
      x = ms_step(x, digits);
      p = p + 1;
    end
    checks(end+1) = value_check('period_from_seed', p, ref, 'period_from_seed');

    bad = double(params.bad_two_digit_targets);
    best = 0;
    best_seed = -1;
    for s = 0:(double(params.two_digit_seeds) - 1)
      x = s;
      n = 0;
      seen2 = false(1, 100);
      while ~any(bad == x) && ~seen2(x + 1)
        seen2(x + 1) = true;
        x = ms_step(x, 2);
        n = n + 1;
      end
      if n > best
        best = n;
        best_seed = s;
      end
    end
    checks(end+1) = value_check('two_digit_max_length', best, ref, ...
                                'two_digit_max_length');
    checks(end+1) = value_check('two_digit_max_length_seed', best_seed, ref, ...
                                'two_digit_max_length_seed');

    res = struct('name', e.name, 'family', e.family, ...
                 'pass', all([checks.pass]), 'checks', checks);
  elseif strcmp(e.name, 'wichmann_hill')
    seeds = double(params.seeds);
    mults = double(params.multipliers);
    mods = double(params.moduli);
    n_samples = double(params.samples);

    periods = zeros(1, 3);
    for i = 1:3
      x = 1;
      seen = false(1, mods(i));
      p = 0;
      while ~seen(x + 1)
        seen(x + 1) = true;
        x = mod(mults(i) * x, mods(i));
        p = p + 1;
      end
      periods(i) = p;
      checks(end+1) = value_check(sprintf('sub_period_%d', i), p, ref, ...
                                  sprintf('sub_period_%d', i));
    end

    lcm = periods(1);
    for i = 2:3
      lcm = lcm * periods(i) / gcd(lcm, periods(i));
    end
    checks(end+1) = value_check('combined_period_lcm', lcm, ref, ...
                                'combined_period_lcm');

    s1 = seeds(1); s2 = seeds(2); s3 = seeds(3);
    p1 = seeds(1); p2 = seeds(2); p3 = seeds(3);
    rs = zeros(1, n_samples);
    max_diff = 0.0;
    for i = 1:n_samples
      s1 = mod(mults(1) * s1, mods(1));
      s2 = mod(mults(2) * s2, mods(2));
      s3 = mod(mults(3) * s3, mods(3));
      r = mod(s1 / mods(1) + s2 / mods(2) + s3 / mods(3), 1);
      rs(i) = r;
      p1 = mults(1) * mod(p1, 177) - 2 * floor(p1 / 177);
      p2 = mults(2) * mod(p2, 176) - 35 * floor(p2 / 176);
      p3 = mults(3) * mod(p3, 178) - 63 * floor(p3 / 178);
      r16 = mod(p1 / mods(1) + p2 / mods(2) + p3 / mods(3), 1);
      max_diff = max(max_diff, abs(r - r16));
      if i <= 4
        checks(end+1) = value_check(sprintf('r_%d', i - 1), r, ref, ...
                                    sprintf('r_%d', i - 1));
      end
    end

    checks(end+1) = value_check('r_min', min(rs), ref, 'r_min');
    checks(end+1) = value_check('r_max', max(rs), ref, 'r_max');
    viol = max(max(0, rs - 1) + max(0, -rs));
    checks(end+1) = value_check('r_range_violation', viol, ref, ...
                                'r_range_violation');
    checks(end+1) = value_check('variant16_max_diff', max_diff, ref, ...
                                'variant16_max_diff');

    res = struct('name', e.name, 'family', e.family, ...
                 'pass', all([checks.pass]), 'checks', checks, ...
                 'combined_period_lcm', lcm, 'r_min', min(rs), ...
                 'r_max', max(rs));
  else
    res = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', ['unknown prng_legacy entry: ' e.name]);
  end
end

function c = value_check(label, computed, ref, refname)
  vref = double(ref.(refname).value);
  tol = double(ref.(refname).abs_tol);
  c = mkcheck(label, computed, vref, tol, ...
              abs(computed - vref) / max(abs(vref), 1), ...
              abs(computed - vref) <= tol);
end

function y = ms_step(x, digits)
  s = sprintf('%0*d', 2 * digits, x * x);
  y = str2double(s((digits / 2 + 1):(digits / 2 + digits)));
end
