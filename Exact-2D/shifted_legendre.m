function y = shifted_legendre(n,a,b,x)
p = legendre(n);
y = 0;
k = n;
for i = 1:n+1
    y = y + p(i) * (2*((x-a)/(b-a))-1).^k;
    k = k-1;
end