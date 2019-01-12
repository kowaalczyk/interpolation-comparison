function c=Bsplnat(y)
  # calculates natural cubic B-spline coefficients
  # input: y - values of a function to interpolate in equally spaced nodes
  # output: c - B-spline coefficients
  n = max(size(y));
  B=build_sparse(n);
  y_comp = padded(y);
  c = B\y_comp;  # if this takes too long, implement Kincaid p. 330 manually
endfunction

function B=build_sparse(n)
  D = sparse(1:n,2:n+1,ones(1,n)*2/3,n,n+2);
  E1 = sparse(1:n,1:n,ones(1,n)/6,n,n+2);
  E2 = sparse(1:n,3:n+2,ones(1,n)/6,n,n+2);
  B = sparse(n+2,n+2);
  B(2:end-1,:) = E1+D+E2;
  B(1,1:3) = [1,-2,1];
  B(end,end-2:end) = [1,-2,1];
endfunction

function yp=padded(y)
  n = max(size(y));
  yp = zeros(n+2,1);
  yp(2:end-1) = y;
endfunction
