# Comparison of interpolation methods

### interpNewton
A function returning Newton-basis coefficients of a polynomial interpolating given points.

### Horner
A vectorized implementation of a Horner algorithm for the polynomials in the NEwton basis.

### Bsplnat
A function calculating natural cubic B-spline coefficients for a given values (must be taken from equidistant points).

### Bsplval
A function evaluating B-spline of a given coefficients in a given vector of points using optimized de Boor algorithm.

### test
A script that demonstrates usege of all functions and tests accuracy of a given methods.
It uses two functions:
- f(x) = cos(3x) for x in [-2pi, pi]
- g(x) = |x| for x in [-10, 10]

to compare:
- polynomial interpolation on equidistant nodes
- polynomial interpolation on Chebyshev nodes (optimal for polynomial interpolation)
- B-spline interpolaion (by definition on equidistant nodes)

using 4, 16 or 64 interpolation nodes.

Resulting interpolation errors are outputted to the console, results are plotted for visual inspection.
