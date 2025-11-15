function [mean, var] = Exact_Moments(f,a,b,c,d,N)

mean = Gauss_Quadrature2D(f,a,b,c,d,N);
g = @ (x,y) (f(x,y)-mean).^2;
var = Gauss_Quadrature2D(g,a,b,c,d,N);