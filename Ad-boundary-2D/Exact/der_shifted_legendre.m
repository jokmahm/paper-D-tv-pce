function dy = der_shifted_legendre(n,a,b,x)
P = legendre(n);
dP = polyder(P);
y = 0;
k = n-1;
for i = 1:n
    y = y + dP(i) * (2*((x-a)/(b-a))-1).^k;
    k = k-1;
end
dy = 2*y;    % chain differentiation