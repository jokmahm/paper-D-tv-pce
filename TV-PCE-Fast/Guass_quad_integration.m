function [Bp, Bp_x, w, V] = Guass_quad_integration(p,a,b,k,f)
% calculate the basis values for 2 D integration

% points and weigths
[P,w] = grule(k);    

% Coordinate transformation
xp(1:k) = 0.5*((b-a).*P(1:k)+(b+a));

% Vectorizedz
Bp = basis(p,a,b,xp);
Bp_x = basis_x(p,a,b,xp);

V = f(xp);