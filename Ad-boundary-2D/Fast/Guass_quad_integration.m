function [B_p, B_x, B_y, W, V] = Guass_quad_integration(p,a,b,c,d,k,f)
% calculate the basis values for 2 D integration

[P,w] = grule(k);    % points and weigths

% Coordinate transformation
xp(1:k) = 0.5*((b-a).*P(1:k)+(b+a));
yp(1:k) = 0.5*((d-c).*P(1:k)+(d+c));

% Vectorized
[Xp,Yp] = meshgrid(xp,yp);
W = w'*w;
B_x = basis_x(p,a,b,c,d,Xp,Yp);
B_y = basis_y(p,a,b,c,d,Xp,Yp);

B_p = basis(p,a,b,c,d,Xp,Yp);
V = f(Xp,Yp);