function out = registry_run(selector)
  % Driver: evaluate registry entries and report computed vs reference.
  % selector: 'all' (default), 'attractors', 'potentials', a family
  % ('map','ifs','ode','neural','potential'), or a single entry name.
  if nargin < 1 || isempty(selector), selector = 'all'; end
  here = fileparts(mfilename('fullpath'));
  root = fileparts(here);
  addpath(here);
  outdir = fullfile(root, 'results');
  if ~exist(outdir, 'dir'), mkdir(outdir); end
  att = registry_load(fullfile(root, 'registries', 'attractors.json'));
  pot = registry_load(fullfile(root, 'registries', 'potentials.json'));
  eng = registry_load(fullfile(root, 'registries', 'energy_drift.json'));
  dsp = registry_load(fullfile(root, 'registries', 'goertzel_fft.json'));
  arc = registry_load(fullfile(root, 'registries', 'archimedes_pi.json'));
  msw = registry_load(fullfile(root, 'registries', 'middle_square_weyl.json'));
  hyr = registry_load(fullfile(root, 'registries', 'hybrid_prng.json'));
  cor = registry_load(fullfile(root, 'registries', 'cordic.json'));
  hast = registry_load(fullfile(root, 'registries', 'hastings.json'));
  lcg = registry_load(fullfile(root, 'registries', 'lcg_lehmer.json'));
  anc = registry_load(fullfile(root, 'registries', 'matematica_ancestral.json'));
  prl = registry_load(fullfile(root, 'registries', 'prng_legacy.json'));

  runs = {};
  for i = 1:numel(att.entries)
    e = att.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      switch e.family
        case 'map',   r = registry_eval_map(e, outdir);
        case 'ifs',   r = registry_eval_ifs(e);
        case 'ode',   r = registry_eval_ode(e);
        case 'neural', r = registry_eval_neural(e);
        otherwise, error('registry_run:family', 'Unknown family %s', e.family);
      end
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end
  for i = 1:numel(pot.entries)
    e = pot.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = registry_eval_potential(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end
  for i = 1:numel(eng.entries)
    e = eng.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = registry_eval_energy(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end

  for i = 1:numel(dsp.entries)
    e = dsp.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = registry_eval_dsp(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end

  for i = 1:numel(arc.entries)
    e = arc.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = registry_eval_archimedes(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end

  for i = 1:numel(msw.entries)
    e = msw.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = registry_eval_middle_square_weyl(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end

  for i = 1:numel(hyr.entries)
    e = hyr.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = registry_eval_hybrid_prng(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end

  runs = run_registry_file(cor, selector, @registry_eval_cordic, runs);
  runs = run_registry_file(hast, selector, @registry_eval_hastings, runs);
  runs = run_registry_file(lcg, selector, @registry_eval_lcg, runs);
  runs = run_registry_file(anc, selector, @registry_eval_ancestral, runs);
  runs = run_registry_file(prl, selector, @registry_eval_prng_legacy, runs);

  npass = 0;
  for i = 1:numel(runs)
    if runs{i}.pass, npass = npass + 1; end
  end
  fprintf('\nSummary: %d/%d entries PASS\n', npass, numel(runs));

  out = struct('timestamp', datestr(now, 30), 'selector', selector, ...
               'n_pass', npass, 'n_total', numel(runs), 'runs', {runs});
  jpath = fullfile(outdir, sprintf('registry_run_%s.json', out.timestamp));
  fid = fopen(jpath, 'w');
  fprintf(fid, '%s', jsonencode(out));
  fclose(fid);
  fprintf('Saved: %s\n', jpath);
end

function runs = run_registry_file(reg, selector, evalfh, runs)
  % Evalua las entradas de un registry con el evaluador indicado.
  for i = 1:numel(reg.entries)
    e = reg.entries{i};
    if ~selected(e, selector), continue; end
    fprintf('[%s] %s\n', e.family, e.name);
    try
      r = evalfh(e);
    catch err
      r = struct('name', e.name, 'family', e.family, 'pass', false, ...
                 'checks', mkcheck('error', NaN, NaN, NaN, NaN, false), ...
                 'error', err.message);
    end
    printrun(r);
    runs{end + 1} = r;
  end
end

function tf = selected(e, sel)
  if strcmp(sel, 'all')
    tf = true;
  elseif strcmp(sel, 'attractors')
    tf = ~strcmp(e.family, 'potential');
  elseif strcmp(sel, 'potentials')
    tf = strcmp(e.family, 'potential');
  else
    tf = strcmp(sel, e.name) || strcmp(sel, e.family);
  end
end

function printrun(r)
  if isfield(r, 'error')
    fprintf('  ERROR: %s\n  => FAIL\n', r.error);
    return;
  end
  for j = 1:numel(r.checks)
    c = r.checks(j);
    if c.pass, tag = 'PASS'; else tag = 'FAIL'; end
    if isnan(c.reference)
      fprintf('  %s %s (observed) = %.12g\n', tag, c.quantity, c.computed);
    else
      fprintf('  %s %s: computed=%.12g ref=%.12g rel_err=%.3g tol=%.3g\n', ...
              tag, c.quantity, c.computed, c.reference, c.rel_err, c.tol);
    end
  end
  if r.pass
    fprintf('  => PASS\n');
  else
    fprintf('  => FAIL\n');
  end
end
