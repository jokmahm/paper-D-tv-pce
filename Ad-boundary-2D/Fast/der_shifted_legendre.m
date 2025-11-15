function dy = der_shifted_legendre(n,a,b,x)

h = 1e-5;
dy = (shifted_legendre(n,a,b,x+h)-shifted_legendre(n,a,b,x))/h;