function phi = basis_y(p,l1,u1,l2,u2,x,y)
global idx
Np = nchoosek(p+2,2);
[m, n] = size(x);
f = zeros(m,n,p+1);
g = zeros(m,n,p+1);
phi = zeros(m,n,Np);

for i = 0:p
    f(:,:,i+1) =  shifted_legendre(i,l1,u1,x);
    g(:,:,i+1) =  der_shifted_legendre(i,l2,u2,y);
end

for i = 1:Np
    phi(:,:,i) =  f(:,:,idx(i,1)+1).*g(:,:,idx(i,2)+1);
end