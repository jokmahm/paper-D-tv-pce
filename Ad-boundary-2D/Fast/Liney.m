function [By_0y, By_1y, By_x0, By_x1, B0y, B1y, w] = Liney(p,a,b,c,d,k)
global idx

[P,w] = grule(k);    % points and weigths
Np = nchoosek(p+2,2);

xp(1:k) = 0.5*((b-a).*P(1:k)+(b+a));
yp(1:k) = 0.5*((d-c).*P(1:k)+(d+c));

B0y = zeros(Np,k);
B1y = zeros(Np,k);

for i = 1:Np
    B0y(i,:) = Phi(i-1,a,b,c,d,a,yp);
    B1y(i,:) = Phi(i-1,a,b,c,d,b,yp);
end

f0 = zeros(k,p+1);
f1 = zeros(k,p+1);
g = zeros(k,p+1);
By_0y = zeros(Np,k);
By_1y = zeros(Np,k);

for i = 0:p
    f0(:,i+1) =  shifted_legendre(i,a,b,a);
    f1(:,i+1) =  shifted_legendre(i,a,b,b);
    g(:,i+1) =  der_shifted_legendre(i,c,d,yp);
end

for i = 1:Np
    By_0y(i,:) =  f0(:,idx(i,1)+1).*g(:,idx(i,2)+1);
    By_1y(i,:) =  f1(:,idx(i,1)+1).*g(:,idx(i,2)+1);
end

f0 = zeros(k,p+1);
f1 = zeros(k,p+1);
g = zeros(k,p+1);
By_x0 = zeros(Np,k);
By_x1 = zeros(Np,k);

for i = 0:p
    f0(:,i+1) =  der_shifted_legendre(i,a,b,c);
    f1(:,i+1) =  der_shifted_legendre(i,a,b,d);
    g(:,i+1) =  shifted_legendre(i,c,d,xp);
end

for i = 1:Np
    By_x0(i,:) =  g(:,idx(i,1)+1).*f0(:,idx(i,2)+1);
    By_x1(i,:) =  g(:,idx(i,1)+1).*f1(:,idx(i,2)+1);
end