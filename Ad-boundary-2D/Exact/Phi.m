function phi = Phi(n,l1,u1,l2,u2,x,y)
global idx
phi = shifted_legendre(idx(n+1,1),l1,u1,x).*shifted_legendre(idx(n+1,2),l2,u2,y);