function v=Horner(c, x, z)
  # vectorized Horner algorithm for the Newton basis
  # input: 
  # c, x are respectively coefficients and nodes of a polynomial in the Newton basis
  # z is the vector of points in which we want to calculate values of the polynomial
  # output: v - values of the polynomial in z
  n = max(size(c));
  v = ones(size(z)) * c(n);
  for j=n-1:-1:1
    v = v .* (z .- x(j)) .+ c(j);
  endfor
endfunction
