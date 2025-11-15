function y = der_shifted_legendre(n,a,b,z)
P = legendre(n);
m = length(P);
dP = (m-1:-1:1).*P(1:m-1);
% dP = polyder(P);
y = 0;
k = n-1;
for i = 1:n
    y = y + dP(i) * (2*((z-a)/(b-a))-1).^k;
    k = k-1;
end
y = 2 * y;   % chain differentiation

% Symbolic form
% p = legendre(n);
% syms x
% g = poly2sym(p);
% g = subs(g,x,2*((x-a)/(b-a))-1);
% dy = diff(g);
% y = subs(dy,x,z);