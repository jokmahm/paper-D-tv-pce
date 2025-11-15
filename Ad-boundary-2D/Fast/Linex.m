function [Bx_x0, Bx_x1, Bx_0y, Bx_1y, Bx0, Bx1, w] = Linex(p,a,b,c,d,k)
global idx

[P,w] = grule(k);    % points and weigths
Np = nchoosek(p+2,2);

xp(1:k) = 0.5*((b-a).*P(1:k)+(b+a));
yp(1:k) = 0.5*((d-c).*P(1:k)+(d+c));

Bx0 = zeros(Np,k);
Bx1 = zeros(Np,k);

for i = 1:Np
    Bx0(i,:) = Phi(i-1,a,b,c,d,xp,c);
    Bx1(i,:) = Phi(i-1,a,b,c,d,xp,d);
end


f0 = zeros(k,p+1);
f1 = zeros(k,p+1);
g = zeros(k,p+1);
Bx_x0 = zeros(Np,k);
Bx_x1 = zeros(Np,k);

for i = 0:p
    f0(:,i+1) =  shifted_legendre(i,c,d,a);
    f1(:,i+1) =  shifted_legendre(i,c,d,b);
    g(:,i+1) =  der_shifted_legendre(i,a,b,xp);
end

for i = 1:Np
    Bx_x0(i,:) =  g(:,idx(i,1)+1).*f0(:,idx(i,2)+1);
    Bx_x1(i,:) =  g(:,idx(i,1)+1).*f1(:,idx(i,2)+1);
end

f0 = zeros(k,p+1);
f1 = zeros(k,p+1);
g = zeros(k,p+1);
Bx_0y = zeros(Np,k);
Bx_1y = zeros(Np,k);

for i = 0:p
    f0(:,i+1) =  der_shifted_legendre(i,c,d,a);
    f1(:,i+1) =  der_shifted_legendre(i,c,d,b);
    g(:,i+1) =   shifted_legendre(i,a,b,yp);
end

for i = 1:Np
    Bx_0y(i,:) =  f0(:,idx(i,1)+1).*g(:,idx(i,2)+1);
    Bx_1y(i,:) =  f1(:,idx(i,1)+1).*g(:,idx(i,2)+1);
end