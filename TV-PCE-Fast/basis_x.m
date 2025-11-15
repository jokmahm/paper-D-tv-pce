function phi = basis_x(p,l,u,x)
n = length(x);
phi = zeros(n,p+1);
for i = 0:p
    phi(:,i+1) =  shifted_legendre_der(i,l,u,x);
end