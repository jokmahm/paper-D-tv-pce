function [mean, var] = Exact_Moments(f,a,b,N)

mean = Gauss_Quadrature1D(f,a,b,N);
g = @ (x) (f(x)-mean).^2;
var = Gauss_Quadrature1D(g,a,b,N);