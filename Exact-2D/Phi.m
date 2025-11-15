function phi = Phi(n,l,u,x,y)
global idx
phi = shifted_legendre(idx(n,1),l,u,x).*shifted_legendre(idx(n,2),l,u,y);