function c=interpNewton(x,y)
  # Based on: Kincaid p. 314
  # input: x, y: f(x) = y where f is the function we wish to interpolate
  # output: Newton basis coefficients for the interpolating polynomial
  c = y;
  n = max(size(x));
  for j = 1:n
    for i = n:-1:j+1
      c(i) = (c(i) - c(i-1)) / (x(i) - x(i-j));
    endfor
  endfor
endfunction
