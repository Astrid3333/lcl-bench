function checks = emptychecks()
  % Empty struct array with the mkcheck field layout, so checks(end+1)= works.
  checks = struct('quantity', {}, 'computed', {}, 'reference', {}, ...
                  'tol', {}, 'rel_err', {}, 'pass', {});
end
