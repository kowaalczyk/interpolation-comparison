# This script demonstrates usage of B-spline and polynomial interpolation functions
# Bsplval, Bsplnat, interpNewton, Horner
# by plotting their results (assignment 3a) and printing interpolation errors (assignment 3b)
# of 2 functions: 
# f(x) = cos(3x) where x in [-2pi, pi]
# g(x) = |x| where x in [-10, 10]
# using varying amount of interpolation nodes

# not a function file:
1;
clear();

N_PLOT_NODES = 1000;  # also used for calculating interpolation error

function y=f(x)
  y = cos(3*x);
endfunction

function y=g(x)
  y = abs(x);
endfunction

function [y_int, y_real, x_plot]=interp_poly_linspace(domain, f, n_nodes, n_pts_to_plot)
  # polynomial interpolation on equidistant nodes prepared for plotting / error calculation
  # input:
  # domain - domain of a function as a pair [domain_start, domain_end]
  # f - function that is to be interpolated
  # n_nodes - number of nodes that will be used for fitting polynomial
  # n_pts_to_plot - number of points for which function values will be generated for plotting
  # output:
  # y_int - interpolating polynomial values calculated in plot points
  # y_real - real function values calculated in plot points
  # x_plot - n_pts_to_plot equidistant plot points
  x = linspace(domain(1), domain(end), n_nodes);
  y = f(x);
  c = interpNewton(x,y);
  x_plot = linspace(domain(1), domain(end), n_pts_to_plot);
  y_real = f(x_plot);
  y_int = Horner(c,x, x_plot);
endfunction

function nodes=Chebyshev(deg, domain=[-1,1])
  # generates roots of Chebyshev's polynomial of degree deg in domain interval
  # based on: https://en.wikipedia.org/wiki/Chebyshev_nodes
  a=domain(1);
  b=domain(end);
  nodes = sort((a+b)/2 + (b-a)/2*cos(pi/deg*((1:deg)-.5)));
endfunction

function [y_int, y_real, x_plot]=interp_poly_chebyshev(domain, f, n_nodes, n_pts_to_plot)
  # polynomial interpolation on Chebyshev nodes prepared for plotting / error calculation
  # input:
  # domain - domain of a function as a pair [domain_start, domain_end]
  # f - function that is to be interpolated
  # n_nodes - number of nodes that will be used for fitting polynomial (=degree of Chebyshev polynomial)
  # n_pts_to_plot - number of points for which function values will be generated for plotting
  # output:
  # y_int - interpolating polynomial values calculated in the plot points
  # y_real - real function values calculated in plot points
  # x_plot - n_pts_to_plot equidistant plot points
  x = Chebyshev(n_nodes, domain);
  y = f(x);
  c = interpNewton(x,y);
  x_plot = linspace(domain(1), domain(end), n_pts_to_plot);
  y_real = f(x_plot);
  y_int = Horner(c,x, x_plot);
endfunction

function [y_int, y_real, x_plot]=interp_bspl(domain, f, n_nodes, n_pts_to_plot)
  # B-spline interpolation on equidistant nodes prepared for plotting / error calculation
  # input:
  # domain - domain of a function as a pair [domain_start, domain_end]
  # f - function that is to be interpolated
  # n_nodes - number of equidistant nodes on which interpolation will be based
  # n_pts_to_plot - number of points for which function values will be generated for plotting
  # output:
  # y_int - interpolating B-spline values calculated in the plot points
  # y_real - real function values calculated in plot points
  # x_plot - n_pts_to_plot equidistant plot points
  x = linspace(domain(1), domain(end), n_nodes);
  y = f(x);
  c = Bsplnat(y);
  x_plot = linspace(domain(1), domain(end), n_pts_to_plot);
  y_real = f(x_plot)';
  y_int = Bsplval(x_plot,c,domain(1), domain(end));
endfunction


# -- assignment 3a ---


# linspace for f
figure();
[y_int1, y_real, x_plot] = interp_poly_linspace([-2*pi(), pi()], @f, 4, N_PLOT_NODES);
y_int2 = interp_poly_linspace([-2*pi(), pi()], @f, 16, N_PLOT_NODES);
plot(x_plot, y_real, x_plot, y_int1, x_plot, y_int2);
title("Polynomial interpolation of f");
legend("true values", "4 nodes linspace", "16 nodes linspace");

# chebyshev for f
figure();
[y_int1, y_real, x_plot] = interp_poly_chebyshev([-2*pi(), pi()], @f, 4, N_PLOT_NODES);
y_int2 = interp_poly_chebyshev([-2*pi(), pi()], @f, 16, N_PLOT_NODES);
plot(x_plot, y_real, x_plot, y_int1, x_plot, y_int2);
title("Polynomial interpolation of f");
legend("true values", "4 nodes Chebyshev", "16 nodes chebyshev");

# bspline for f
figure();
[y_int1, y_real, x_plot] = interp_bspl([-2*pi(), pi()], @f, 4, N_PLOT_NODES);
y_int2 = interp_bspl([-2*pi(), pi()], @f, 16, N_PLOT_NODES);
plot(x_plot, y_real, x_plot, y_int1, x_plot, y_int2);
title("B-spline interpolation of f");
legend("true values", "4 nodes linspace", "16 nodes linspace");

# linspace for g
figure();
[y_int1, y_real, x_plot] = interp_poly_linspace([-10, 10], @g, 4, N_PLOT_NODES);
y_int2 = interp_poly_linspace([-10, 10], @g, 16, N_PLOT_NODES);
plot(x_plot, y_real, x_plot, y_int1, x_plot, y_int2);
title("Polynomial interpolation of g");
legend("true values", "4 nodes linspace", "16 nodes linspace");

# chebyshev for g
figure();
[y_int1, y_real, x_plot] = interp_poly_chebyshev([-10, 10], @g, 4, N_PLOT_NODES);
y_int2 = interp_poly_chebyshev([-10,10], @g, 16, N_PLOT_NODES);
plot(x_plot, y_real, x_plot, y_int1, x_plot, y_int2);
title("Polynomial interpolation of g");
legend("true values", "4 nodes Chebyshev", "16 nodes Chebyshev");

# bspline for g
figure();
[y_int1, y_real, x_plot] = interp_bspl([-10, 10], @g, 4, N_PLOT_NODES);
y_int2 = interp_bspl([-10,10], @g, 16, N_PLOT_NODES);
plot(x_plot, y_real, x_plot, y_int1, x_plot, y_int2);
title("B-spline interpolation of g");
legend("true values", "4 nodes linspace", "16 nodes linspace");


# -- assignment 3b ---


for n=[4,16,64]
  printf("Interpolation erors using %d nodes:\n", fix(n));
  # linspace for f
  [y_int, y_real, x_plot] = interp_poly_linspace([-2*pi(), pi()], @f, n, N_PLOT_NODES);
  printf("Linspace for f: %g", norm(y_int - y_real, Inf));
  printf("\n");
  # chebyshev for f
  [y_int, y_real, x_plot] = interp_poly_chebyshev([-2*pi(), pi()], @f, n, N_PLOT_NODES);
  printf("Chebyshev for f: %g", norm(y_int - y_real, Inf));
  printf("\n");
  # bspline for f
  [y_int, y_real, x_plot] = interp_bspl([-2*pi(), pi()], @f, n, N_PLOT_NODES);
  printf("Bspline for f: %g", norm(y_int - y_real, Inf));
  printf("\n");
  # linspace for g
  [y_int, y_real, x_plot] = interp_poly_linspace([-10, 10], @g, n, N_PLOT_NODES);
  printf("Linspace for g: %g", norm(y_int - y_real, Inf));
  printf("\n");
  # chebyshev for g
  [y_int, y_real, x_plot] = interp_poly_chebyshev([-10, 10], @g, n, N_PLOT_NODES);
  printf("Chebyshev for g: %g", norm(y_int - y_real, Inf));
  printf("\n");
  # bspline for g
  [y_int, y_real, x_plot] = interp_bspl([-10, 10], @g, n, N_PLOT_NODES);
  printf("Bspline for g: %g", norm(y_int - y_real, Inf));
  printf("\n\n");
endfor
