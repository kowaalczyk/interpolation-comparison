function v = Bsplval(z, c, a, b)
  # calculate B-spline values
  # input:
  # z - vector of points to calculate values on
  # c - B-spline coefficients (created using Bsplnat)
  # a, b - boundaries of an interval on which B-spline equidistant nodes are placed
  # output:
  # v - vector of B-spline values at z
  t = build_padded_knot_pos(a, b, max(size(c))-2);
  n_nodes = max(size(z));
  v = zeros(n_nodes, 1);
  for i = 1 : n_nodes
    v(i) = de_boor(z(i),t,c);
  endfor
endfunction

function t=build_padded_knot_pos(a,b,n_knots)
  raw_t = linspace(a, b, n_knots)';
  dist = raw_t(2)-raw_t(1);  # distance between any 2 consecutive knots
  t = [linspace(a-3*dist, a-dist, 3)'; raw_t; linspace(b+dist, b+3*dist, 3)'];
endfunction

function v=de_boor(x,padded_knot_pos,control_points)
  # calculate value for a single point x, using optimized de Boore algorithm,
  # based on: https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
  knot_interval_idx = find_interval(x, padded_knot_pos);
  d = control_points(knot_interval_idx-3:knot_interval_idx);
  for r = 1:3
    for j = 3:-1:r
      alpha = (x - padded_knot_pos(j+knot_interval_idx-3)) ...
          / (padded_knot_pos(j+1+knot_interval_idx-r) ...
          - padded_knot_pos(j+knot_interval_idx-3));
      d(j+1) = (1-alpha)*d(j) + alpha*d(j + 1);
    endfor
  endfor
  v = d(end);
endfunction

function k=find_interval(z,tab)
  # finds smallest i that satisfies tab[i] <= z < tab[i+1]
  # if such i does not exist (z is too small or too large),
  # returns the corresponding border of the unpadded tab
  n = length(tab);
  k = 0;
  for i = 1 : (n - 1)
    k = i;
    if tab(i) <= z && tab(i+1) > z
      break;
    endif
  endfor
  # make sure not to exceed the padding (normalize i to the border of unpadded tab):
  k = min(k,n-4);
  k = max(k,3);
endfunction
